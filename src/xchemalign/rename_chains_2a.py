# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# This is a hacky one-off script that handles renaming of some 2A protease structures where some structures are
# homo-dimers, but the chain naming is inconsistent - some use chains A and B and some use A and C.
# This script renames those that are A and C to become A and B, with chain B (usually the Zinc atom) being renamed
# to chain C.
# This script is specific for this case and not for general use.


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
        xtal = row.get(utils.Constants.SOAKDB_XTAL_NAME)
        pdb = row.get(utils.Constants.SOAKDB_COL_PDB)
        cif = row.get(utils.Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST)
        pdb_handled = False
        cif_handled = False
        if pdb:
            before, key, after = pdb.partition('model_building/')
            if after:
                pdb_p = mb_p / after
                if pdb_p.is_file():
                    print(index, pdb_p)
                    struc = gemmi.read_pdb(str(pdb_p), ignore_ter=True)
                    struc.setup_entities()
                    rename_chains(xtal, pdb_p, struc)
                # else:
                #     print(pdb_p, 'not found')
                pdb_handled = True
            if not pdb_handled:
                print('FAILED TO HANDLE', index, pdb)

        if cif:
            before, key, after = cif.partition('model_building/')
            if after:
                cif_p = mb_p / after
                if cif_p.is_file():
                    print(index, cif_p)
                    struc = gemmi.read_structure(str(cif_p))
                    struc.setup_entities()
                    rename_chains(xtal, cif_p, struc)
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
        rename_chains(ref, Path(ref), struc)


def rename_chains(xtal, pth, struc):
    model = struc[0]
    chain_b = None
    chain_c = None
    res_b = 0
    res_c = 0
    pol_a = 0
    pol_b = 0
    pol_c = 0
    for chain in model:
        if chain.name == 'A':
            pol_a = chain.get_polymer().length()
            print('polymer A:', chain.get_polymer().length())
        elif chain.name == 'B':
            chain_b = chain
            res_b = len(chain)
            pol_b = chain.get_polymer().length()
            print('polymer B:', chain.get_polymer().length())
        elif chain.name == 'C':
            chain_c = chain
            res_c = len(chain)
            pol_c = chain.get_polymer().length()
            print('polymer C:', chain.get_polymer().length())

    # assume that protein chains are long and hetatm ones are short
    if pol_a > 125 and pol_b > 125:
        print(xtal, 'has expected A and B chains')
    elif pol_c > 125 and pol_b > 0:
        print('ERROR: chain B polymer length = {} and chain C polymer length = {}'.format(pol_b, pol_c))
    elif pol_c > 125:
        dump_chains(model)

        pth.rename(pth.parent / (pth.name + '.orig'))

        chain_c.name = 'B'
        if chain_b:
            chain_b.name = 'C'

        # this renames the subchains (asym id) to be consistent with the chain names
        struc.assign_subchains(force=True)

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
