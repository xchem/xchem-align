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


def split_altloc_pdb(pdb_path: Path, ligand_names: list[str]):
    """
    Take a PDB file as input and generate a PDB file for each ligand altloc
    :return:
    """

    # First pass - find out if there are any ligands with altcodes
    with open(pdb_path, 'rt') as pdb:
        ligand_alt_codes = set()
        for line in pdb:
            if line.startswith('HETATM'):
                residue_name = line[17:21].strip()
                if residue_name in ligand_names:
                    alt_code = line[16:17]
                    if alt_code != ' ':
                        ligand_alt_codes.add(alt_code)

        print(ligand_alt_codes)

        for i, ligand_alt_code in enumerate(ligand_alt_codes):
            # second pass - split out the PDB files using those altcodes
            with open(pdb_path, 'rt') as pdb:
                # generate the output file:
                dir = pdb_path.parent
                base_filename = pdb_path.stem
                print(dir, base_filename)
                with open(dir / (base_filename + '_' + ligand_alt_code + '.pdb'), 'wt') as out:
                    for line in pdb:
                        if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('ANISOU'):
                            alt_code = line[16:17]
                            if alt_code in ligand_alt_codes and alt_code != ligand_alt_code:
                                print("SKIPPING ", ligand_alt_code, i, line.strip())
                                continue
                        out.write(line)


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
        self.ligand_name = None

    def validate(self):
        errors = 0
        if not self.pdbfile:
            self.log("PDB file not defined:", level=2)
            errors += 1
        if not os.path.isfile(self.pdbfile):
            self.log("PDB file not found:", self.pdbfile, level=2)
            errors += 1
        if self.biomol and not os.path.isfile(self.biomol):
            self.log("biomol file not found:", self.biomol, level=2)
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
            self.log("attaching biomol headers")
            self.add_biomol_remark()

    def create_apo_solv_desolv(self):
        """
        Creates two files:
        _apo-desolv - as _apo, but without solvent, ions and buffers;
        _apo-solv - just the ions, solvent and buffers

        :returns: Created file names
        """

        if not self.apo_file:
            self.log("apo file has not been created. Use pdb_apo().make_apo_file()")
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

    def _extract_residue_as_list(self, chain: str, res_id: int):
        """
        Extracts out the lines corresponding a particular chain and residue
        :param chain: e.g. A
        :param res_id: e.g. 123
        :return: The lines of the PDB corresponding to the specified residue
        """

        ligand_lines = []
        with open(self.pdbfile, "r") as pdb:
            for line in pdb:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    if line[21] == chain and line[22:26].strip() == res_id:
                        ligand_lines.append(line)
        return ligand_lines

    def extract_residue(self, chain, res_id):
        return ''.join(self._extract_residue_as_list(chain, res_id))

    def extract_coordinates(self, chain, res_id):
        lines = self._extract_residue_as_list(chain, res_id)
        name = lines[0][17:21].strip()
        result = {}
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                atom_id = line[12:16].strip()
                altloc = line[16]
                if altloc == ' ':
                    altloc = '_'
                d = result.get(altloc)
                if d is None:
                    d = {}
                    result[altloc] = d
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                d[atom_id] = (x, y, z)
        return result, name

    def create_ligands(self, chain: str, res_id, cif_file: str):
        """
        Generate a Molfile for the specified ligand
        :param chain: e.g. A
        :param res_id: e.g. 123 (int or str)
        :param cif_file: CIF file with the expected topology
        :return: The generated RDKit Mol
        """

        self.log("creating ligand files")

        # get the coordinates of the ligand corresponding to the specified residue number
        coords_for_altlocs, name = self.extract_coordinates(chain, res_id)
        if not coords_for_altlocs:
            self.log('no coordinates found for residue ' + chain + '/' + str(res_id), level=1)
            return None

        # read the ligands from the CIF
        mols = utils.gen_mols_from_cif(cif_file)

        mol = None
        for m in mols:
            if name == m.GetProp('_Name'):
                mol = m
                break
        if not mol:
            self.log('no molecule named ' + name + ' in CIF file', level=1)
            return None

        mol.RemoveAllConformers()
        generated_mols = []

        old_name = mol.GetProp('_Name')
        self.ligand_name = old_name
        self.ligand_base_file = self.output_dir / (self.filebase + "_ligand_" + old_name)
        self.smiles = Chem.MolToSmiles(mol)

        for altloc, coords in coords_for_altlocs.items():
            m = Chem.RWMol(mol)
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
                self.log('atoms not present in cif molecule: ' + str(missing), level=1)

            m.SetProp('_Name', str(self.filebase) + '_' + altloc + old_name)
            m.AddConformer(conf)
            Chem.AssignStereochemistryFrom3D(m)

            modified_mol = self.handle_covalent_ligands(chain, res_id, m)
            if modified_mol:
                m = modified_mol
            generated_mols.append(m)

        with Chem.SDWriter(str(self.ligand_base_file) + '.sdf') as writer:
            for m in generated_mols:
                writer.write(m)

        merged_mol = Chem.RWMol()
        for m in generated_mols:
            merged_mol.InsertMol(m)

        Chem.MolToMolFile(merged_mol, str(self.ligand_base_file) + '.mol')

        # write the ligand PDB
        with open(str(self.ligand_base_file) + '.pdb', 'wt') as pdb_out:
            with open(self.pdbfile) as pdb_in:
                for line in pdb_in:
                    if line.startswith('HETATM') or line.startswith('ANISOU'):
                        # match residue number, name and chain
                        if (
                            line[22:26].strip() == str(res_id)
                            and line[17:21].strip() == old_name
                            and chain == line[21].strip()
                        ):
                            pdb_out.write(line)

        with open(str(self.ligand_base_file) + '.smi', 'wt') as f:
            f.write(self.smiles)

        return merged_mol

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

    def handle_covalent_ligands(self, chain: str, res_num: str, lig_mol: Chem.Mol):
        """
        Look for LINK records that define a covalent bond between the ligand and the protein
        e.g.
        LINK         SG  CYS A 110                 C7  LIG A 201     1555   1555  2.22

        If found then create a new Mol that has the protein atom added to it and a bond to that atom

        :param chain:
        :return: A modified molecule if one was needed, if not None
        """

        link_records = []
        L_ATOM_TYPE = 'l_atom_type'
        L_RES_NAME = 'l_res_name'
        L_RES_NUM = 'l_res_num'
        L_CHAIN = 'l_chain'
        P_ATOM_TYPE = 'p_atom_type'
        P_RES_NAME = 'p_res_name'
        P_RES_NUM = 'p_res_num'
        P_CHAIN = 'p_chain'

        mol = None

        with open(self.pdbfile, "r") as pdb:
            header_complete = False
            for line in pdb:
                if not header_complete:
                    if line.startswith("ATOM"):
                        header_complete = True
                    elif line.startswith("LINK"):
                        # The ligand could be first or second in the link line. We need to work out which.
                        # We do this by matching the residue ids
                        left_atom_type = line[12:16].strip()  # chars 13-16
                        left_res_name = line[17:20].strip()  # chars 18-20
                        left_chain = line[21]  # char 22
                        left_res_num = line[22:26].strip()  # chars 23-26

                        right_atom_type = line[42:46].strip()  # chars 43-46
                        right_res_name = line[47:50].strip()  # chars 48-50
                        right_chain = line[51]  # char 52
                        right_res_num = line[52:56].strip()  # chars 53-56

                        # now set the ligand and protein properties
                        found_ligand = True
                        d = {}
                        if res_num == left_res_num:
                            d[L_ATOM_TYPE] = left_atom_type
                            d[L_RES_NAME] = left_res_name
                            d[L_CHAIN] = left_chain
                            d[L_RES_NUM] = left_res_num
                            d[P_ATOM_TYPE] = right_atom_type
                            d[P_RES_NAME] = right_res_name
                            d[P_CHAIN] = right_chain
                            d[P_RES_NUM] = right_res_num
                            link_records.append(d)
                        elif res_num == right_res_num:
                            d[L_ATOM_TYPE] = right_atom_type
                            d[L_RES_NAME] = right_res_name
                            d[L_CHAIN] = right_chain
                            d[L_RES_NUM] = right_res_num
                            d[P_ATOM_TYPE] = left_atom_type
                            d[P_RES_NAME] = left_res_name
                            d[P_CHAIN] = left_chain
                            d[P_RES_NUM] = left_res_num
                            link_records.append(d)
                        else:
                            found_ligand = False
                            self.log(
                                "failed to identify ligand residue. Looking for",
                                res_num,
                                "found",
                                left_res_num,
                                right_res_num,
                            )

                        # if found_ligand and chain == d[L_CHAIN]:
                        #     self.log("|".join([d[P_ATOM_TYPE], d[P_RES_NAME], d[P_CHAIN], d[P_RES_NUM],
                        #                        d[L_ATOM_TYPE], d[L_RES_NAME], d[L_CHAIN], d[L_RES_NUM],]))

                else:
                    # look for the protein atom
                    if line.startswith("ATOM"):
                        for link_record in link_records:
                            if (
                                line[23:26].strip() == link_record[P_RES_NUM]
                                and line[21] == link_record[P_CHAIN]
                                and line[12:16].strip() == link_record[P_ATOM_TYPE]
                            ):
                                x = line[30:38].strip()  # chars 31-38
                                y = line[38:46].strip()  # chars 39-46
                                z = line[46:54].strip()  # chars 47-54
                                element = line[76:78].strip()  # chars 77-78

                                if len(element) == 2:
                                    element[1] = element[1].lower()
                                atomic_no = Chem.GetPeriodicTable().GetAtomicNumber(element)

                                self.log(
                                    "found covalent atom",
                                    "|".join(
                                        [
                                            d[P_ATOM_TYPE],
                                            d[P_RES_NAME],
                                            d[P_CHAIN],
                                            d[P_RES_NUM],
                                            d[L_ATOM_TYPE],
                                            d[L_RES_NAME],
                                            d[L_CHAIN],
                                            d[L_RES_NUM],
                                        ]
                                    ),
                                    element,
                                    atomic_no,
                                    x,
                                    y,
                                    z,
                                )
                                mol = Chem.RWMol(lig_mol)
                                # find the ligand atom to attach the bond to
                                lig_atom = None
                                for atom in mol.GetAtoms():
                                    if atom.GetProp('atom_id').upper() == link_record[L_ATOM_TYPE].upper():
                                        lig_atom = atom
                                        break
                                if lig_atom:
                                    idx = mol.AddAtom(Chem.Atom(atomic_no))
                                    mol.AddBond(idx, lig_atom.GetIdx(), Chem.BondType.SINGLE)
                                    mol.GetConformer(-1).SetAtomPosition(
                                        idx, Geometry.Point3D(float(x), float(y), float(z))
                                    )
                                    Chem.SanitizeMol(mol)
                                    Chem.AssignStereochemistryFrom3D(mol)
                                    self.log("added atom and bond for covalent ligand")
                                else:
                                    self.log("could not find ligand atom", link_record[L_ATOM_TYPE], level=1)
                                break

        return mol

    def log(self, *args, level=0):
        if self.logger:
            self.logger.log(*args, level=level)
        else:
            print(*args)


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
