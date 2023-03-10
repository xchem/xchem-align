import argparse, os, json
from pathlib import Path


class PDBXtal:

    def __init__(self, pdbfile, output_dir, biomol=None):
        self.pdbfile = pdbfile
        self.filebase = Path(pdbfile).stem
        self.output_dir = output_dir
        self.biomol = biomol
        self.non_ligs = json.load(open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r"))
        self.apo_file = None
        print('filebase:', self.filebase)

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

    def make_apo_file(self, keep_headers=False):
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

        self.apo_file = os.path.join(self.output_dir, str(self.filebase + "_apo.pdb"))
        print('reading:', self.apo_file)
        pdb = open(self.apo_file, "w")
        pdb.write(str(lines))
        pdb.close()

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

        prot_file_name = os.path.join(self.output_dir, str(self.filebase + "_apo-desolv.pdb"))
        prot_file = open(prot_file_name, "w")
        solv_file_name = os.path.join(self.output_dir, str(self.filebase + "_apo-solv.pdb"))
        solv_file = open(solv_file_name, "w")

        for line in open(self.apo_file).readlines():
            if line.startswith("HETATM"):
                solv_file.write(line)
            else:
                prot_file.write(line)
        solv_file.close()
        prot_file.close()

        return prot_file_name, solv_file_name

    def create_pdb_for_ligand(self, ligand, count, reduce, smiles_file, covalent=False):
        """
        A pdb file is produced for an individual ligand, containing atomic and connection information
        :param ligand: Name of the Ligand
        :param count: The index of the ligand
        :param reduce: Bool, if the file needs to be named using the chain name of the PDB
        :param smiles_file: File path of smiles_file (if any)
        :param covalent: Bool, indicate whether or not covalent attach should be sought.
        :return: .pdb file for ligand.
        """
        output_file_name = os.path.join(self.output_dir, self.filebase + '_ligand.pdb')

        individual_ligand = []
        individual_ligand_conect = []
        # adding atom information for each specific ligand to a list
        for atom in self.final_hets:
            if str(atom[16:20].strip() + atom[20:26]) == str(ligand):
                individual_ligand.append(atom)

        con_num = 0
        for atom in individual_ligand:
            atom_number = atom.split()[1]
            for conection in self.conects:
                if (atom_number in conection and conection not in individual_ligand_conect):
                    individual_ligand_conect.append(conection)
                    con_num += 1

        # checking that the number of conect files and number of atoms are almost the same
        # (taking into account ligands that are covalently bound to the protein

        # assert 0 <= con_num - len(individual_ligand) <= 1
        # making into one list that is compatible with conversion to mol object
        ligand_het_con = individual_ligand + individual_ligand_conect
        # make a pdb file for the ligand molecule

        if os.path.isdir(lig_out_dir):
            # This is stupid but will correctly spec the files... is there a better solution??
            shutil.rmtree(lig_out_dir)

        if not os.path.isdir(lig_out_dir):
            os.makedirs(lig_out_dir)

        ligands_connections = open(
            os.path.join(lig_out_dir, (file_base + ".pdb")), "w+"
        )
        for line in ligand_het_con:
            ligands_connections.write(str(line))
        ligands_connections.close()
        # making pdb file into mol object
        mol = self.create_pdb_mol(
            file_base=file_base, lig_out_dir=lig_out_dir, smiles_file=smiles_file, handle_cov=covalent)
        # Move Map files into lig_out_dir

        if not mol:
            print(
                f'WARNING: {file_base} did not produce a mol object from its pdb lig file!')
        else:
            try:
                Chem.AddHs(mol)

                self.mol_lst.append(mol)
                self.mol_dict["directory"].append(lig_out_dir)
                self.mol_dict["mol"].append(mol)
                self.mol_dict["file_base"].append(file_base)

            except AssertionError:
                print(file_base, 'is unable to produce a ligand file')
                pass

    def create_mol_file(self, mol_obj, smiles_file=None):
        """
        a .mol file is produced for an individual ligand

        :param mol_obj: The RDKit Mol file object
        :param smiles_file: The filepath of a text file that contains the smiles string of the mol file (if exists).
        :return: A mol file!
        """

        out_file = os.path.join(directory, str(file_base + ".mol"))

        if not mol_obj:
            print(f'WARNING: mol object is empty: {file_base}')

        if smiles_file:
            with open(smiles_file, 'r') as sf:
                smiles_list = sf.readlines()
            mol_dicts = {}
            sim_dicts = {}
            original_fp = Chem.RDKFingerprint(mol_obj)
            for smiles in smiles_list:
                try:
                    template = AllChem.MolFromSmiles(smiles.rstrip())
                    new_mol = AllChem.AssignBondOrdersFromTemplate(template, mol_obj)
                    new_fp = Chem.RDKFingerprint(new_mol)
                    mol_dicts[smiles] = new_mol
                    sim_dicts[smiles] = DataStructs.FingerprintSimilarity(original_fp, new_fp)
                except Exception as e:
                    print(e)
                    print('failed to fit template ' + smiles_file)
                    print(f'template smiles: {smiles}')
            # If none of the templates fit...
            # Just write the smiles and use the mol_obj
            # If there is 1, just use it
            # If there are multiple that fit, find the one with the highest similarity
            if len(mol_dicts) < 1:
                mol_obj = mol_obj
            else:
                mol_obj = mol_dicts[max(sim_dicts, key=sim_dicts.get)]
                smiles = max(sim_dicts, key=sim_dicts.get)
        else:
            print(f'Warning: No smiles file: {file_base}')
            smiles = Chem.MolToSmiles(mol_obj)

        # Write output mol file...
        Chem.rdmolfiles.MolToMolFile(mol_obj, out_file)

        # Write new smiles_txt
        smiles_out_file = os.path.join(directory, str(file_base + "_smiles.txt"))
        with open(smiles_out_file, 'w+') as smiles_txt:
            smiles_txt.write(smiles)
        # Create output png too
        out_png = os.path.join(directory, str(file_base + ".png"))
        draw_mol = Chem.Mol(mol_obj)
        AllChem.Compute2DCoords(draw_mol)
        Draw.MolToFile(draw_mol, out_png, imageType='png')
        return mol_obj


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

    p.make_apo_file(keep_headers=args.keep_headers)
    p.create_apo_solv_desolv()
    print('Done')


if __name__ == "__main__":
    main()
