import numpy as np
from loguru import logger

from ligand_neighbourhood_alignment import dt 
from ligand_neighbourhood_alignment.structure import _get_dataset_protein_chains

def _get_closest_xtalform(xtalforms: dict[str, dt.XtalForm], structure, structures):
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    for xtalform_id, xtalform in xtalforms.items():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        # Check they are in the same spacegroup
        if ref_spacegroup != structure_spacegroup:
            continue

        # Check they have the same protein chains
        xtalform_protein_chains = [
            _chain for xtalform_assembly in xtalform.assemblies.values() for _chain in xtalform_assembly.chains
        ]
        dataset_protein_chains = _get_dataset_protein_chains(structure)

        if set(dataset_protein_chains) != set(xtalform_protein_chains):
            continue

        deltas = np.array(
            [
                structure_cell.a / ref_structure_cell.a,
                structure_cell.b / ref_structure_cell.b,
                structure_cell.c / ref_structure_cell.c,
                structure_cell.alpha / ref_structure_cell.alpha,
                structure_cell.beta / ref_structure_cell.beta,
                structure_cell.gamma / ref_structure_cell.gamma,
            ]
        )
        xtalform_deltas[xtalform_id] = deltas

    if len(xtalform_deltas) == 0:
        return None, None

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform]




def _assign_dataset(dataset, assemblies, xtalforms, structure, structures):
    closest_xtalform_id, deltas = _get_closest_xtalform(
        xtalforms,
        structure,
        structures,
    )

    if (closest_xtalform_id is None) & (deltas is None):
        logger.info(f"No reference in same spacegroup as: {dataset.dtag}")
        logger.info(f"Structure path is: {dataset.pdb}")
        raise Exception(f"No reference in same spacegroup as: {dataset.dtag}\nStructure path is: {dataset.pdb}")

    if np.any(deltas > 1.1) | np.any(deltas < 0.9):
        logger.info(f"No reference for dataset: {dataset.dtag}")
        logger.info(f"Deltas to closest unit cell are: {deltas}")
        logger.info(f"Structure path is: {dataset.pdb}")

        raise Exception(
            f"No reference for dataset: {dataset.dtag}\nDeltas to closest unit cell in {closest_xtalform_id} are: {deltas}\nStructure path is: {dataset.pdb}"
        )

    return closest_xtalform_id