import numpy as np
from loguru import logger
import gemmi
from rich import print as print
import itertools

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
        ref_centroids[chain] = get_chain_centroid(ref[0][chain])
        assert ref_centroids[chain].size == 3

    # Get the mov chain centroids
    mov_centroids = {}
    # Get an alignment based on the first chain (should get around indexing choice!)
    ref_poly = ref[0][xtalform_protein_chains[0]].get_polymer()
    mov_poly = mov[0][xtalform_protein_chains[0]].get_polymer()
    sup = gemmi.calculate_superposition(
        ref_poly, 
        mov_poly, 
        ref_poly.check_polymer_type(), 
        gemmi.SupSelect.CaP,
        )
    
    for chain in xtalform_protein_chains:
        transformed_chain = mov[0][chain].clone().get_polymer()
        transformed_chain.transform_pos_and_adp(sup.transform)

        mov_centroids[chain] = get_chain_centroid(transformed_chain)
        assert mov_centroids[chain].size == 3

    # Get the distances under symmetry and PBC
    distances = {}
    # for indexing in itertools.product([-1, 1], [-1, 1], [-1])

    for ref_chain, ref_centroid in ref_centroids.items():
        distances[ref_chain] = {}
        for mov_chain, mov_centroid in mov_centroids.items():
            nearest_image = ref.cell.find_nearest_image(
                gemmi.Position(ref_centroid[0], ref_centroid[1], ref_centroid[2]), 
                gemmi.Position(mov_centroid[0], mov_centroid[1], mov_centroid[2]), 
                gemmi.Asu.Any,
                )
            
            dist = nearest_image.dist()
            print(f'Closest distance between {ref_centroid} and {mov_centroid} is {dist}')
            distances[ref_chain][mov_chain] = dist
            
    # Assign each movin chain its closest ref chain
    assignments = {}
    min_distances = {}
    for ref_chain in xtalform_protein_chains:
        closest_chain = min(
            distances[ref_chain], 
            key=lambda x: distances[ref_chain][x],
            )
        min_distances[ref_chain] = distances[ref_chain][closest_chain]
        assignments[ref_chain] = closest_chain

    if len(set(assignments.values())) != len(set(assignments.keys())):
        print(f'ERROR! Chain assignment is {assignments}. No two chains should be modelled in the same symmetry position in the unit cell!')
        raise Exception()
    
    return assignments, min_distances, ref_centroids, mov_centroids



    ...

def _get_closest_xtalform(xtalforms: dict[str, dt.XtalForm], structure, structures):
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    all_distances = {}
    all_mappings = {}
    all_ref_centroids = {}
    all_mov_centroids = {}
    for xtalform_id, xtalform in xtalforms.items():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        # Check they are in the same spacegroup
        if ref_spacegroup != structure_spacegroup:
            continue

        # Check they have the same protein chains
        xtalform_protein_chains = [
            _chain for xtalform_assembly in xtalform.assemblies.values() 
            for _chain in xtalform_assembly.chains
        ]
        dataset_protein_chains = _get_dataset_protein_chains(structure)

        if set(dataset_protein_chains) != set(xtalform_protein_chains):
            continue

        # Check if the chains align reasonably
        chain_mapping, min_distances, ref_centroids, mov_centroids = get_xtalform_chain_mapping(
            ref_structure, 
            structure, 
            xtalform_protein_chains,
            )
        print(f'Chanin mapping: {chain_mapping}')
        print(f'Min distances: {min_distances}')
        all_mappings[xtalform_id] = chain_mapping
        all_distances[xtalform_id] = min_distances
        all_ref_centroids[xtalform_id] = ref_centroids
        all_mov_centroids[xtalform_id] = mov_centroids

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
        return None, None, all_mappings, all_distances, all_ref_centroids, all_mov_centroids

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform], all_mappings, all_distances, all_ref_centroids, all_mov_centroids




def _assign_dataset(dataset, assemblies, xtalforms, structure, structures):
    closest_xtalform_id, deltas, all_mappings, all_distances, all_ref_centroids, all_mov_centroids = _get_closest_xtalform(
        xtalforms,
        structure,
        structures,
    )

    if (closest_xtalform_id is None) & (deltas is None):
        xtalform_info = {
            xtalform_id: structures[xtalform.reference].spacegroup_hm
            for xtalform_id, xtalform in xtalforms.items()
        }
        message = (
            f"No reference for: {dataset.dtag}\n"
            f"Structure path is: {dataset.pdb}\n"
            'Usually the best solution is to create a new crystalform with this dataset as a reference.\n'
            'However there -are- several possible problems:\n'
            '1. There is a genuine missing crystalform - for example none with the same spacegroup\n'
            '   and/or unit cell as this dataset. This can be diagnosed by comparing the \n'
            '   "# Dataset spacegroup" entry to the "# xtalform spacegroups" entry below. This is fixed\n'
            '   by creating a new xtalform with this dataset as reference.\n'
            '2. Chain names vary to the references in a way the program cannot deal with. \n' 
            '   In particular if after alignment on one chain -other- chains are not in approximately\n' \
            '   the same place as their counterparts in a reference. This can be diagnosed by comparing\n' \
            '   seeing if all chains map to themselves in the "# All chain mappings" entry below. This can\n' \
            '   be fixed by renaming chains in this dataset or removing this dataset.\n'
            '# Diagnostic Information for someone who knows what they are doing:\n'
            '# All chain mappings\n'
            f'{all_mappings}\n'
            '# All Distances\n'
            f'{all_distances}\n'
            '# All ref centroids\n'
            f'{all_ref_centroids}\n'
            '# All mov centroids\n'
            f'{all_mov_centroids}\n'
            f'# Dataset spacegroup: {structure.spacegroup_hm}\n'
            '# Xtalform spacegroups\n'
            f'{xtalform_info}\n'
        )
        raise Exception(message)

    if np.any(deltas > 1.1) | np.any(deltas < 0.9):
        logger.info(f"No reference for dataset: {dataset.dtag}")
        logger.info(f"Deltas to closest unit cell are: {deltas}")
        logger.info(f"Structure path is: {dataset.pdb}")

        raise Exception(
            f"No reference for dataset: {dataset.dtag}\nDeltas to closest unit cell in {closest_xtalform_id} are: {deltas}\nStructure path is: {dataset.pdb}"
        )

    return closest_xtalform_id