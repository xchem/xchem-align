import argparse

from pathlib import Path

import gemmi

from xchemalign import dbreader
from xchemalign import utils


def run(visit_dir, polymer_length=100):
    visit_p = Path(visit_dir)
    db_p = visit_p / 'processing/database/soakDBDataFile.sqlite'
    mb_p = visit_p / 'processing/analysis/model_building'

    name, df = dbreader.read_all(db_p)

    print('SoakDB dataframe dimensions:', name, df.shape)

    for index, row in df.iterrows():
        xtal = row.get(utils.Constants.SOAKDB_XTAL_NAME)
        pdb = row.get(utils.Constants.SOAKDB_COL_PDB)
        handled = False
        if pdb:
            before, key, after = pdb.partition('model_building/')
            if after:
                pdb_p = mb_p / after
                if pdb_p.is_file():
                    struc = gemmi.read_pdb(str(pdb_p), ignore_ter=True)
                    struc.setup_entities()
                    list_chains(xtal, struc, polymer_length)
                handled = True
    if pdb and not handled:
        print('FAILED TO HANDLE', index, pdb)


def list_chains(xtal, struc, polymer_length):
    model = struc[0]
    chains = ''
    for chain in model:
        if chain.get_polymer().length() >= polymer_length:
            chains += chain.name.lower()

    if chains:
        print('        - ' + xtal + '  ' + chains)


def main():
    # example:
    # python -m xchemalign.list_chains -d data/std_test/lb32627-66_2026-02-13/inputs/dls/labxchem/data/lb32627/lb32627-66
    #
    # To report structures that only have polymers of length 75 or greater in chains A and B use something like this:
    # python -m xchemalign.list_chains -d data/std_test/lb32627-66_2026-02-13/inputs/dls/labxchem/data/lb32627/lb32627-66 -p 75 | grep a | grep b | cut -c 1-23
    # This outputs text that can be C&P'd into the config.yaml
    # Note: chain names are output in lower case so that they can be grepped (assuming the crystal names are always upper case)

    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-d", "--dir", help="visit dir")
    parser.add_argument("-p", "--polymer-length", type=int, default=100, help="min length of polymer")
    args = parser.parse_args()

    run(args.dir, polymer_length=args.polymer_length)


if __name__ == "__main__":
    main()
