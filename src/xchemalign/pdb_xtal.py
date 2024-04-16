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

from rdkit import Chem
from rdkit import Geometry

from xchemalign import utils


class PDBXtal:
    def __init__(self, pdbfile, output_dir, biomol=None, logger=None):
        self.pdbfile = pdbfile
        self.filebase = Path(pdbfile).stem
        self.output_dir = Path(output_dir)
        self.biomol = biomol
        self.logger = logger
        self.non_ligs = json.load(open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r"))
        self.apo_file = None
        self.ligand_file = None
        self.apo_solv_file = None
        self.apo_desolv_file = None
        self.ligand_base_file = None
        self.smiles = None

    def validate(self):
        errors = 0
        if not self.pdbfile:
            self.log("PDB file not defined:", level=2)
            errors += 1
        if not os.path.isfile(self.pdbfile):
            self.log("PDB file not found:", self.pdbfile, level=2)
            errors += 1
        if self.biomol and not os.path.isfile(self.biomol):
            self.log("Biomol file not found:", self.biomol, level=2)
            errors += 1

        return errors

    def add_biomol_remark(self):
        """
        Add contents of biomol/additional text file to the _apo.pdb file.
        This results in re-writing the _apo.pdb file.
        """
        biomol_remark = open(self.biomol).read().rstrip() + "\n"

        f = self.apo_file
        with open(f) as handle:
            switch = 0
            header_front, header_end = [], []
            pdb = []
            for line in handle:
                if line.startswith("ATOM"):
                    switch = 1
                if line.startswith("HETATM"):
                    switch = 2
                if switch == 0:
                    header_front.append(line)
                elif (switch == 2) and not line.startswith("HETATM"):
                    header_end.append(line)
                else:
                    pdb.append(line)
            full_file = "".join(header_front) + biomol_remark + "".join(pdb) + "".join(header_end)
            with open(f, "w") as w:
                w.write(full_file)

    def create_apo_file(self, keep_headers=False):
        """
        Keeps anything other than unique ligands

        :param: pdb file
        :returns: created XXX_apo.pdb file
        """
        lines = ""
        ligand = ""
        conect_lines = ""
        if keep_headers:
            include = ["CONECT", "SEQRES", "TITLE", "ANISOU"]
        else:
            include = ["CONECT", "REMARK", "CRYST", "SEQRES", "HEADER", "TITLE", "ANISOU"]

        with open(self.pdbfile, "r") as pdb:
            found_atoms = False
            num_lig_lines = 0

            for line in pdb:
                if found_atoms:
                    if line.startswith("CONECT"):
                        conect_lines += line
                    if line[17:20] == "LIG":
                        num_lig_lines += 1
                        ligand += line
                    else:
                        lines += line
                else:
                    if line.startswith("ATOM"):
                        found_atoms = True
                    lines += line

        # prune the conect records to only include those from the LIG
        lig_atom_nums = set()
        conect_to_keep = ""
        for line in ligand:
            tokens = line.split()
            if line.startswith("HETATM"):
                lig_atom_nums.add(int(tokens[1]))
        for line in conect_lines:
            tokens = line.split()
            missing = False
            for token in tokens[1:]:
                if int(token) not in lig_atom_nums:
                    missing = True
                    continue
            if missing:
                print("discarding", line)
            else:
                conect_to_keep += line

        self.apo_file = self.output_dir / (self.filebase + "_apo.pdb")
        self.ligand_file = self.output_dir / (self.filebase + "_ligand_0.pdb")
        self.ligand_file = self.output_dir / (self.filebase + "_ligand_0.pdb")
        with open(self.apo_file, "w") as f:
            f.write(str(lines))
        with open(self.ligand_file, "w") as f:
            f.write(str(ligand))
            f.write(str(conect_lines))

        if self.biomol is not None:
            self.log("Attaching biomol headers")
            self.add_biomol_remark()

    def create_apo_solv_desolv(self):
        """
        Creates two files:
        _apo-desolv - as _apo, but without solvent, ions and buffers;
        _apo-solv - just the ions, solvent and buffers

        :returns: Created file names
        """

        if not self.apo_file:
            self.log("Apo file has not been created. Use pdb_apo().make_apo_file()")
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

    def _extract_residue_as_list(self, chain, res_id):
        """
        Extracts out the lines corresponding a particular chain and residue
        :param chain: e.g. A
        :param res_id: e.g. 123
        :return: The lines of the PDB corresponding to the specified residue
        """

        s = str(res_id)
        if len(s) == 1:
            id = '   ' + s
        elif len(s) == 2:
            id = '  ' + s
        elif len(s) == 3:
            id = ' ' + s
        elif len(s) == 4:
            id = s
        else:
            raise ValueError('Unexpected residue ID: ' + res_id)

        ligand_lines = []
        with open(self.pdbfile, "r") as pdb:
            for line in pdb:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    if line[21] == chain and line[22:26] == id:
                        ligand_lines.append(line)
        return ligand_lines

    def extract_residue(self, chain, res_id):
        return ''.join(self._extract_residue_as_list(chain, res_id))

    def extract_coordinates(self, chain, res_id):
        lines = self._extract_residue_as_list(chain, res_id)
        result = {}
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                atom_id = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                result[atom_id] = (x, y, z)
        return result

    def create_ligands(self, chain: str, res_id, cif_file: str):
        """
        Generate a Molfile for the specified ligand
        :param chain: e.g. A
        :param res_id: e.g. 123 (int or str)
        :param cif_file: CIF file with the expected topology
        :return: The generated RDKit Mol
        """

        mol = utils.gen_mol_from_cif(cif_file)
        mol.RemoveAllConformers()

        coords = self.extract_coordinates(chain, res_id)
        conf = Chem.Conformer()
        missing = []
        for atom_id, pos in coords.items():
            found = False
            for atom in mol.GetAtoms():
                if atom.GetProp('atom_id').upper() == atom_id.upper():
                    idx = atom.GetIdx()
                    point = Geometry.Point3D(pos[0], pos[1], pos[2])
                    conf.SetAtomPosition(idx, point)
                    found = True
                    continue
            if not found:
                missing.append(atom_id)
        if missing:
            self.log('atoms not present in cif molecule:', missing, level=1)

        mol.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(mol)

        self.ligand_base_file = self.output_dir / (self.filebase + "_ligand")
        mol.SetProp('_Name', str(self.filebase))
        self.smiles = Chem.MolToSmiles(mol)
        Chem.MolToMolFile(mol, str(self.ligand_base_file) + '.mol')
        Chem.MolToPDBFile(mol, str(self.ligand_base_file) + '.pdb')
        with open(str(self.ligand_base_file) + '.smi', 'wt') as f:
            f.write(self.smiles)

        return mol

    def extract_sequences(self, pdb_file=None):
        if not pdb_file:
            pdb_file = self.pdbfile

        with open(pdb_file, "rt") as pdb:
            lines = pdb.readlines()
            curr_chain = None
            curr_resno = 0
            curr_seq = None
            seqs = []
            for line in lines:
                if line.startswith('ATOM'):
                    alt = line[16].strip()
                    chain = line[21].strip()
                    code = line[17:20].strip()
                    resno = int(line[22:26].strip())
                    if chain != curr_chain:
                        curr_chain = chain
                        curr_seq = ProteinSeq(chain, [], start=int(resno))
                        seqs.append(curr_seq)
                    if resno != curr_resno:
                        for i in range(resno - curr_resno - 1):
                            curr_seq.seq.append('UNK')
                        curr_resno = resno
                        curr_seq.seq.append(code)
            return seqs

    def log(self, *args, level=0):
        if self.logger:
            self.logger.log(args, level=level)
        else:
            print(args)


class ProteinSeq:
    def __init__(self, chain, seq, start=1):
        self.chain = chain
        self.seq = seq
        self.start = start

    def create_seqres_header(self):
        line_res = []
        lines = []
        line_count = 1
        num_res_str = str(len(self.seq)).rjust(4, ' ')
        for i, res in enumerate(self.seq):
            if i > 0 and i % 13 == 0:
                lines.append(
                    'SEQRES {} {} {}  {}'.format(
                        str(line_count).rjust(3, ' '), self.chain, num_res_str, ' '.join(line_res)
                    )
                )
                line_count += 1
                line_res = []
            line_res.append(res)
        if line_res:
            lines.append(
                'SEQRES {} {} {}  {}'.format(
                    str(line_count).rjust(3, ' '), self.chain, num_res_str, ' '.join(line_res)
                )
            )
        return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description="pdb_xtal")

    parser.add_argument("-p", "--pdb", required=True, help="PDB file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-k", "--keep-headers", action="store_true", help="Keep headers")
    parser.add_argument("-b", "--biomol_txt", help="Biomol Input txt file")
    parser.add_argument("-c", "--chain", default="A", help="Protein chain")
    parser.add_argument("-r", "--residue", help="Ligand residue number")
    parser.add_argument("-f", "--cif", help="CIF file with the ligand topology")

    args = parser.parse_args()
    print("pdb_xtal: ", args)

    p = PDBXtal(args.pdb, args.output_dir, biomol=args.biomol_txt)
    errors = p.validate()
    if errors:
        print("There were validation errors. Please fix and re-run")
        exit(1)

    p.create_apo_file(keep_headers=args.keep_headers)
    p.create_apo_solv_desolv()
    if args.residue and args.cif:
        lig_mol = p.create_ligands(args.chain, args.residue, args.cif)
        print(lig_mol)
    # seqs = p.extract_sequences()
    # for seq in seqs:
    #     print(seq.create_seqres_header())
    print("Done")


if __name__ == "__main__":
    main()
