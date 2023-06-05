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

from rdkit import Chem

from standardize_molecule import standardize_to_iso_smiles


def molecule_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)


def standardize_molecule(smi):
    """
    Generate the standard SMILES for this molecule

    :param smi: The SMILES to standardize
    :return: The standardized SMILES and the RDKit Mol
    """
    std_smi, mol = standardize_to_iso_smiles(smi)
    return std_smi, mol


def generate_inchi(mol, opts=''):
    inchis = Chem.inchi.MolToInchi(mol, options=opts)
    inchik = Chem.inchi.InchiToInchiKey(inchis)
    return inchis, inchik


def remove_stereochemistry(mol):
    mol.RemoveStereochemistry()


def main():

    smi = 'N[C@@H](C)C(=O)O'
    std_smi, mol = standardize_molecule(smi)
    inchis1, inchik1 = generate_inchi(mol)
    Chem.rdmolops.RemoveStereochemistry(mol)
    inchis2, inchik2 = generate_inchi(mol)
    print(smi, std_smi, inchis1, inchik1, inchis2, inchik2)


if __name__ == "__main__":
    main()