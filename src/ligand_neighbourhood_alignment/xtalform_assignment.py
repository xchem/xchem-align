import numpy as np
from loguru import logger
import gemmi

from ligand_neighbourhood_alignment import dt 
from ligand_neighbourhood_alignment.structure import _get_dataset_protein_chains

def chain_to_array(chain):
    poss = []
    for res in chain:
        for atom in res:
            if "CA" not in atom.name:
                continue
            pos = atom.pos
            poss.append(
                [
                    pos.x,
                    pos.y,
                    pos.z
                ]
            )

    return np.array(poss)
    ...

def get_chain_centroid(chain):
    chain_array = chain_to_array(chain)
    print(chain_array.shape)
    return np.mean(chain_array, axis=0)
    ...

def get_xtalform_chain_mapping(ref, mov, xtalform_protein_chains):
    """
    For ref each chain, get the centroid
    For each moving chain, get the centroid
    Get the symmetric distance from each moving chain to each reference chain
    Get the asignment with the minimum rmsd
    """

    # Get the ref chain centroids
    ref_centroids = {}
    for chain in xtalform_protein_chains:
        get_chain_centroid(ref[0][chain])

    # Get the mov chain centroids
    mov_centroids = {}
    for chain in xtalform_protein_chains:
        get_chain_centroid(mov[0][chain])

    # Get the distances under symmetry and PBC
    distances = {}
    for ref_chain, ref_centroid in ref_centroids.items():
        distances[ref_chain] = {}
        for mov_chain, mov_centroid in mov_centroids.items():
            nearest_image = ref.cell.find_nearest_image(
                gemmi.Position(ref_centroid[0], ref_centroid[1], ref_centroid[2]), 
                gemmi.Position(mov_centroid[0], mov_centroid[1], mov_centroid[2]), 
                gemmi.Asu.Any,
                )
            dist = nearest_image.dist()
            distances[ref_chain][mov_chain] = dist
            
    # Assign each movin chain its closest ref chain
    assignments = {}
    min_distances = {}
    for ref_chain in assignments:
        closest_chain = min(assignments[ref_chain], key=lambda x: assignments[ref_chain][x])
        min_distances[ref_chain] = distances[ref_chain][closest_chain]
        assignments[ref_chain] = closest_chain

    if len(set(assignments.values())) != len(set(assignments.keys())):
        print(f'ERROR! Chain assignment is {assignments}. No two chains should be modelled in the same symmetry position in the unit cell!')
        raise Exception()
    
    return assignments, min_distances



    ...

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

        # Check if the chains align reasonably
        chain_mapping, min_distances = get_xtalform_chain_mapping(ref_structure, structure, xtalform_protein_chains)
        print(chain_mapping)
        print(min_distances)

        if not all([x == chain_mapping[x] for x in chain_mapping]):
            print('Chain Mapping is degenerate!')
            print(chain_mapping)
            continue

        if np.mean([x for x in min_distances.values()]) > 5:
            print('Min distances is large!')
            print(min_distances)
            continue

        # Get the deltas
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