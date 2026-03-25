import argparse
import glob
from pathlib import Path
import gemmi

from xchemalign import dbreader
from xchemalign import utils


# print(dir(gemmi))


def dump_chains(model):
    for chain in model:
        print(chain, chain.name)


def run(visit_dir):
    visit_p = Path(visit_dir)
    db_p = visit_p / 'processing/database/soakDBDataFile.sqlite'
    mb_p = visit_p / 'processing/analysis/model_building'
    ref_p = visit_p / 'processing/reference'

    name, df = dbreader.read_all(db_p)

    print('DF:', name, df.shape)

    for index, row in df.iterrows():
        pdb = row.get(utils.Constants.SOAKDB_COL_PDB)
        cif = row.get(utils.Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST)
        pdb_handled = False
        cif_handled = False
        if pdb:
            before, key, after = pdb.partition('model_building/')
            if after:
                pdb_p = mb_p / after
                print(index, pdb_p)
                if pdb_p.is_file():
                    struc = gemmi.read_pdb(str(pdb_p), ignore_ter=True)
                    struc.setup_entities()
                    rename_chains(pdb_p, struc)
                # else:
                #     print(pdb_p, 'not found')
                pdb_handled = True
            if not pdb_handled:
                print('FAILED TO HANDLE', index, pdb)

        if cif:
            before, key, after = cif.partition('model_building/')
            if after:
                cif_p = mb_p / after
                print(index, cif_p)
                if cif_p.is_file():
                    struc = gemmi.read_structure(str(cif_p))
                    struc.setup_entities()
                    rename_chains(cif_p, struc)
                # else:
                #     print(cif_p, 'not found')
                cif_handled = True
            if not cif_handled:
                print('FAILED TO HANDLE', index, cif)

    refs = glob.glob(str(ref_p) + '/*.pdb')
    for ref in refs:
        print('handling ref', ref)
        struc = gemmi.read_pdb(ref, ignore_ter=True)
        struc.setup_entities()
        rename_chains(Path(ref), struc)


def rename_chains(pth, struc):
    model = struc[0]
    chain_b = None
    chain_c = None
    res_b = 0
    res_c = 0
    for chain in model:
        if chain.name == 'B':
            chain_b = chain
            res_b = len(chain)
        if chain.name == 'C':
            chain_c = chain
            res_c = len(chain)

    # assume that protein chains are long and hetatm ones are short
    if res_b < 100 and res_c > 100:
        dump_chains(model)
        pth.rename(pth.parent / (pth.name + '.orig'))

        chain_c.name = 'B'
        if chain_b:
            chain_b.name = 'C'

        if str(pth).endswith('.pdb'):
            struc.write_pdb(str(pth))
        elif str(pth).endswith('.cif'):
            doc = struc.make_mmcif_document()
            doc.write_file(str(pth))


def main():
    # example:
    # python -m xchemalign.rename_chains_2a -d data/std_test/lb32627-66_2026-02-13/inputs/dls/labxchem/data/lb32627/lb32627-66
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-d", "--dir", help="visit dir")
    args = parser.parse_args()

    run(args.dir)


if __name__ == "__main__":
    main()
