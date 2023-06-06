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

import argparse
import json
import os
from pathlib import Path


class PDBXtal:

    def __init__(self, pdbfile, output_dir, biomol=None):
        self.pdbfile = pdbfile
        self.filebase = Path(pdbfile).stem
        self.output_dir = output_dir
        self.biomol = biomol
        self.non_ligs = json.load(open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r"))
        self.apo_file = None
        self.apo_solv_file = None
        self.apo_desolv_file = None

    def validate(self):
        errors = 0
        if not self.pdbfile:
            print('PDB file not defined:')
            errors += 1
        if not os.path.isfile(self.pdbfile):
            print('PDB file not found:', self.pdbfile)
            errors += 1
        if self.biomol and not os.path.isfile(self.biomol):
            print('Biomol file not found:', self.biomol)
            errors += 1

        return errors

    def add_biomol_remark(self):
        """
        Add contents of biomol/additional text file to the _apo.pdb file.
        This results in re-writing the _apo.pdb file.
        """
        biomol_remark = open(self.biomol).read().rstrip() + '\n'

        f = self.apo_file
        with open(f) as handle:
            switch = 0
            header_front, header_end = [], []
            pdb = []
            for line in handle:
                if line.startswith('ATOM'):
                    switch = 1
                if line.startswith('HETATM'):
                    switch = 2
                if switch == 0:
                    header_front.append(line)
                elif (switch == 2) and not line.startswith('HETATM'):
                    header_end.append(line)
                else:
                    pdb.append(line)
            full_file = ''.join(header_front) + biomol_remark + ''.join(pdb) + ''.join(header_end)
            with open(f, 'w') as w:
                w.write(full_file)

    def create_apo_file(self, keep_headers=False):
        """
        Keeps anything other than unique ligands

        :param: pdb file
        :returns: created XXX_apo.pdb file
        """
        lines = ""
        if keep_headers:
            include = ['CONECT', 'SEQRES', 'TITLE', 'ANISOU']
        else:
            include = ['CONECT', 'REMARK', 'CRYST', 'SEQRES', 'HEADER', 'TITLE', 'ANISOU']

        with open(self.pdbfile, 'r') as pdb:
            for line in pdb:
                if (
                        line.startswith("HETATM")
                        and line.split()[3] not in self.non_ligs
                        or any([line.startswith(x) for x in include])
                ):
                    continue
                else:
                    lines += line

        self.apo_file = self.output_dir / (self.filebase + "_apo.pdb")
        f = open(self.apo_file, "w")
        f.write(str(lines))
        f.close()

        if self.biomol is not None:
            print('Attaching biomol headers')
            self.add_biomol_remark()

    def create_apo_solv_desolv(self):
        """
        Creates two files:
        _apo-desolv - as _apo, but without solvent, ions and buffers;
        _apo-solv - just the ions, solvent and buffers

        :returns: Created file names
        """

        if not self.apo_file:
            print( "Apo file has not been created. Use pdb_apo().make_apo_file()")
            exit(1)

        prot_file_name = self.output_dir / (self.filebase + "_apo-desolv.pdb")
        prot_file = open(prot_file_name, "w")
        self.apo_desolv_file = prot_file_name
        solv_file_name = self.output_dir / (self.filebase + "_apo-solv.pdb")
        solv_file = open(solv_file_name, "w")
        self.apo_solv_file = solv_file_name

        for line in open(self.apo_file).readlines():
            if line.startswith("HETATM"):
                solv_file.write(line)
            else:
                prot_file.write(line)
        solv_file.close()
        prot_file.close()

        return prot_file_name, solv_file_name


def main():

    parser = argparse.ArgumentParser(description='pdb_xtal')

    parser.add_argument('-p', '--pdb', required=True, help='PDB file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-k', "--keep-headers", action='store_true', help='Keep headers')
    parser.add_argument('-b', '--biomol_txt', help='Biomol Input txt file')

    args = parser.parse_args()
    print("pdb_xtal: ", args)

    p = PDBXtal(args.pdb, args.output_dir, biomol=args.biomol_txt)
    errors = p.validate()
    if errors:
        print('There were validation errors. Please fix and re-run')
        exit(1)

    p.create_apo_file(keep_headers=args.keep_headers)
    p.create_apo_solv_desolv()
    print('Done')


if __name__ == "__main__":
    main()
