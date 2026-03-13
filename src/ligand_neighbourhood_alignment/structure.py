from loguru import logger
import gemmi

from ligand_neighbourhood_alignment import constants

def _get_dataset_protein_chains(structure):
    print('Check if ZN in res names')
    print('ZN' in constants.RESIDUE_NAMES)
    protein_chains = []
    for model in structure:
        for chain in model:
            if any([(_residue.name in constants.RESIDUE_NAMES) for _residue in chain]):
                protein_chains.append(chain.name)

    return protein_chains

def get_res_from_resid(resid, structure):
    return structure[0][resid[0]][resid[1]][0]
