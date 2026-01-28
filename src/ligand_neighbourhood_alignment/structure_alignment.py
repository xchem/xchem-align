# import os
import re
from pathlib import Path

import gemmi
import networkx as nx
from loguru import logger
import numpy as np
from rich import print as rprint

from ligand_neighbourhood_alignment import dt


from ligand_neighbourhood_alignment.dt import (
    AssignedXtalForms,
    CanonicalSites,
    ConformerSites,
    LigandNeighbourhoods,
    Output,
    SiteTransforms,
    Transforms,
    XtalForms,
    gemmi_to_transform,
    transform_to_gemmi,
)
from ligand_neighbourhood_alignment import alignment_heirarchy
from ligand_neighbourhood_alignment.alignment_landmarks import icode_to_string


def superpose_structure(transform, structure):
    new_structure = structure.clone()

    for model in new_structure:
        for chain in model:
            span = chain.whole()
            span.transform_pos_and_adp(transform)

    return new_structure


def align_structure(
    _structure,
    moving_ligand_id,
    reference_ligand_id,
    g,
    transforms,
    site_transforms: SiteTransforms,
    canonical_site_id,
    conformer_site_id,
    out_path,
):
    shortest_path = nx.shortest_path(g, moving_ligand_id, reference_ligand_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = moving_ligand_id
    running_transform = gemmi.Transform()
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = transforms.get_transform(
                (
                    next_ligand_id,
                    previous_ligand_id,
                ),
            )
            running_transform = transform.combine(running_transform)

        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Subsite alignment transform
    subsite_transform = transform_to_gemmi(
        site_transforms.get_conformer_site_transform(canonical_site_id, conformer_site_id)
    )
    subsite_transform.combine(running_transform)

    # Site alignment transform
    site_transform = transform_to_gemmi(site_transforms.get_canonical_site_transform(canonical_site_id))
    site_transform.combine(running_transform)

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(running_transform)}")

    _structure = superpose_structure(running_transform, _structure)

    # Write the fully aligned structure
    _structure.write_pdb(str(out_path))


def _mark_atom_pos_to_ni_pos_tup(point, mark, st):
    cra = mark.to_cra(st[0])
    ni = st.cell.find_nearest_pbc_position(point, cra.atom.pos, mark.image_idx)
    return (ni.x, ni.y, ni.z)


def _drop_non_binding_chains_and_symmetrize_waters(
    _structure, neighbourhood, moving_ligand_id, dataset_ligand_neighbourhood_ids, xtalform, xtalform_sites
):
    # Other Ligand IDs
    other_ligand_ids = [
        (lid[1], lid[2])
        for lid in dataset_ligand_neighbourhood_ids
        if not ((lid[1] == moving_ligand_id[1]) & (lid[2] == moving_ligand_id[2]))
    ]
    # Get a copy of structure to edit
    new_structure = _structure.clone()

    # Determine which chains have binding residues
    neighbourhood_chains = set([_atom_id[0] for _atom_id in neighbourhood.atoms])

    # Determine the assembly each chain is part of
    chain_assemblies = {
        _chain: _assembly for _assembly_name, _assembly in xtalform.assemblies.items() for _chain in _assembly.chains
    }

    # Get the assembly the ligand is modelled as part of
    for xsid, _xtalform_site in xtalform_sites.items():
        _xtalform_id = _xtalform_site.xtalform_id
        if moving_ligand_id in _xtalform_site.members:
            xtalform_site = _xtalform_site
    site_chain = xtalform_site.crystallographic_chain
    lig_assembly = chain_assemblies[site_chain]

    # Determine which waters are bound near the ligand, and at what positions
    ns = gemmi.NeighborSearch(new_structure[0], new_structure.cell, 10).populate(include_h=False)

    # # Iterate ligand heavy atoms, finding marks and collating
    all_marks = {}
    atom_multiplicities = {}
    for atom in new_structure[0][moving_ligand_id[1]][moving_ligand_id[2]][0]:
        point = atom.pos
        marks = ns.find_atoms(point, '\0', radius=5)

        # # Get cras
        for mark in marks:
            cra = mark.to_cra(new_structure[0])

            if cra.residue.name != 'HOH':
                continue

            base_atom_id = (
                cra.chain.name,
                str(cra.residue.seqid.num)+icode_to_string(cra.residue.seqid.icode),
                cra.atom.name,
            )

            # Note first occurence of each atom
            if base_atom_id not in atom_multiplicities:
                atom_multiplicities[base_atom_id] = 1

            # Get the atom id with multiplicity
            atom_id = base_atom_id + (atom_multiplicities[base_atom_id],)

            pos1 = _mark_atom_pos_to_ni_pos_tup(point, mark, new_structure)

            #
            if atom_id in all_marks:
                # Get current marks with same base_atom_id
                comparator_ids = [(a, b, c, d) for a, b, c, d in all_marks if (a, b, c) == base_atom_id]

                # Check if it is distinct from all of them
                poss = [
                    _mark_atom_pos_to_ni_pos_tup(
                        all_marks[comparator_id][0], all_marks[comparator_id][1], new_structure
                    )
                    for comparator_id in comparator_ids
                ]

                # If so increase the multiplicity by one and add the atom with the new multiplicity
                if not any([np.allclose(pos1, pos2, atol=0.1) for pos2 in poss]):
                    atom_multiplicities[base_atom_id] += 1
                    atom_id = base_atom_id + (atom_multiplicities[base_atom_id],)

                # Otherwise skip
                else:
                    continue

            all_marks[atom_id] = (point, mark, pos1)

    # Update water positions if they are near ligand but modelled elsewhere
    local_water_chains = {}
    chain_name_to_chain = {_chain.name: _chain for _chain in new_structure[0]}
    for atom_id, (point, mark, mark_pos) in all_marks.items():
        # Get the corresponding atom
        cra = mark.to_cra(new_structure[0])

        # if not a water, skip
        if cra.residue.name != 'HOH':
            continue

        # If a symatom, find the symchain to associate it with
        if atom_id[3] > 1:
            chain_name = f'{cra.chain.name}{atom_id[3]}'
            # # If associated with a new symchain, add it to the structure
            if chain_name not in chain_name_to_chain:
                chain = gemmi.Chain(chain_name)
                new_structure[0].add_chain(chain)
                chain_name_to_chain[chain_name] = chain

            # # Otherwise get the knewn chain
            chain = new_structure[0][chain_name]

            # Add the new residue, and select the relevant atom
            residue = cra.residue.clone()
            chain.add_residue(residue)
            residue = new_structure[0][chain_name][atom_id[1]][0]
            atom = residue[atom_id[2]][0]

        # Otherwise get the original atom for modification
        else:
            chain = cra.chain
            residue = cra.residue
            atom = cra.atom

        # If water update position and note chain
        if residue.name == 'HOH':
            # Record local water chain and seqid
            if chain.name not in local_water_chains:
                local_water_chains[chain.name] = []
            local_water_chains[chain.name].append(str(residue.seqid.num)+icode_to_string(residue.seqid.icode))

            # Update water position from mark
            pos = atom.pos
            pos.x = mark_pos[0]
            pos.y = mark_pos[1]
            pos.z = mark_pos[2]

    # Drop residues and non-local waters from non-binding chains containing site waters
    new_chains = []
    for _model in new_structure:
        for _chain in _model:
            # Create a new chain to hold symmetry waters from non-binding chains
            new_chain = gemmi.Chain(_chain.name)
            # Iterate residues in the old chain, adding the local waters
            for _residue in _chain:
                if (_chain.name, str(_residue.seqid.num)+icode_to_string(_residue.seqid.icode)) in other_ligand_ids:
                    continue
                if _residue.name == 'HOH':
                    if _chain.name in local_water_chains:
                        if str(_residue.seqid.num)+icode_to_string(residue.seqid.icode) in local_water_chains[_chain.name]:
                            new_chain.add_residue(_residue.clone())
                else:
                    if (_chain.name in lig_assembly.chains) or (_chain.name in neighbourhood_chains):
                        # Don't include other ligands
                        new_chain.add_residue(_residue.clone())

            if len(new_chain) > 0:
                new_chains.append(new_chain)

    for new_chain in new_chains:
        del new_structure[0][new_chain.name]
        new_structure[0].add_chain(new_chain)

    return new_structure


def _drop_non_assembly_chains_and_symmetrize_waters(
    _structure, 
    neighbourhood, 
    moving_ligand_id, 
    dataset_ligand_neighbourhood_ids, xtalform, xtalform_sites
):
    DEBUG = False
    if (moving_ligand_id[0] in ['6w7h', '6w8q']) :
        print(f'Processing error dataset!')
        DEBUG = True

    # Other Ligand IDs
    other_ligand_ids = [
        (lid[1], lid[2])
        for lid in dataset_ligand_neighbourhood_ids
        if not ((lid[1] == moving_ligand_id[1]) & (lid[2] == moving_ligand_id[2]) )
    ]
    if DEBUG:
        print(f'OTHER LIGAND IDS: {moving_ligand_id} : {other_ligand_ids}')
    # Get a copy of structure to edit
    new_structure = _structure.clone()

    # Determine which chains have binding residues
    neighbourhood_chains = set([_atom_id[0] for _atom_id in neighbourhood.atoms])
    if DEBUG:
        print(f'NEIGHBOURHOOD CHAINS: {moving_ligand_id} : {neighbourhood_chains}')

    # Determine the assembly each chain is part of
    chain_assemblies = {
        _chain: _assembly for _assembly_name, _assembly in xtalform.assemblies.items() for _chain in _assembly.chains
    }
    if DEBUG:
        print(f'CHAIN ASSEMBLY: {moving_ligand_id} : {chain_assemblies}')
                    

    # Get the assembly the ligand is modelled as part of
    for xsid, _xtalform_site in xtalform_sites.items():
        _xtalform_id = _xtalform_site.xtalform_id
        if moving_ligand_id in _xtalform_site.members:
            xtalform_site = _xtalform_site
    site_chain = xtalform_site.crystallographic_chain
    lig_assembly = chain_assemblies[site_chain]
    if DEBUG:
        print(f'SITE CHAIN: {moving_ligand_id} : {site_chain}')
        print(f'LIGAND ASSEMBLY: {moving_ligand_id} : {lig_assembly}')

    # Determine which waters are bound near the ligand, and at what positions
    ns = gemmi.NeighborSearch(new_structure[0], new_structure.cell, 10).populate(include_h=False)

    # # Iterate ligand heavy atoms, finding marks and collating
    all_marks = {}
    atom_multiplicities = {}
    for atom in new_structure[0][moving_ligand_id[1]][moving_ligand_id[2]][0]:
        point = atom.pos
        marks = ns.find_atoms(point, '\0', radius=5)

        # # Get cras
        for mark in marks:
            cra = mark.to_cra(new_structure[0])

            if cra.residue.name != 'HOH':
                continue

            base_atom_id = (
                cra.chain.name,
                str(cra.residue.seqid.num)+icode_to_string(cra.residue.seqid.icode),
                cra.atom.name,
            )

            # Note first occurence of each atom
            if base_atom_id not in atom_multiplicities:
                atom_multiplicities[base_atom_id] = 1

            # Get the atom id with multiplicity
            atom_id = base_atom_id + (atom_multiplicities[base_atom_id],)

            pos1 = _mark_atom_pos_to_ni_pos_tup(point, mark, new_structure)

            #
            if atom_id in all_marks:
                # Get current marks with same base_atom_id
                comparator_ids = [(a, b, c, d) for a, b, c, d in all_marks if (a, b, c) == base_atom_id]

                # Check if it is distinct from all of them
                poss = [
                    _mark_atom_pos_to_ni_pos_tup(
                        all_marks[comparator_id][0], all_marks[comparator_id][1], new_structure
                    )
                    for comparator_id in comparator_ids
                ]

                # If so increase the multiplicity by one and add the atom with the new multiplicity
                if not any([np.allclose(pos1, pos2, atol=0.1) for pos2 in poss]):
                    atom_multiplicities[base_atom_id] += 1
                    atom_id = base_atom_id + (atom_multiplicities[base_atom_id],)

                # Otherwise skip
                else:
                    continue

            all_marks[atom_id] = (point, mark, pos1)

    # Update water positions if they are near ligand but modelled elsewhere
    local_water_chains = {}
    chain_name_to_chain = {_chain.name: _chain for _chain in new_structure[0]}

    for atom_id, (point, mark, mark_pos) in all_marks.items():
        # Get the corresponding atom
        cra = mark.to_cra(new_structure[0])

        # if not a water, skip
        if cra.residue.name != 'HOH':
            continue

        # If a symatom, find the symchain to associate it with
        if atom_id[3] > 1:
            chain_name = f'{cra.chain.name}{atom_id[3]}'
            # # If associated with a new symchain, add it to the structure
            if chain_name not in chain_name_to_chain:
                chain = gemmi.Chain(chain_name)
                new_structure[0].add_chain(chain)
                chain_name_to_chain[chain_name] = chain

            # # Otherwise get the knewn chain
            chain = new_structure[0][chain_name]

            # Add the new residue, and select the relevant atom
            residue = cra.residue.clone()
            chain.add_residue(residue)
            residue = new_structure[0][chain_name][atom_id[1]][0]
            atom = residue[atom_id[2]][0]

        # Otherwise get the original atom for modification
        else:
            chain = cra.chain
            residue = cra.residue
            atom = cra.atom

        # If water update position and note chain
        if residue.name == 'HOH':
            # Record local water chain and seqid
            if chain.name not in local_water_chains:
                local_water_chains[chain.name] = []
            local_water_chains[chain.name].append(str(residue.seqid.num)+icode_to_string(residue.seqid.icode))

            # Update water position from mark
            pos = atom.pos
            pos.x = mark_pos[0]
            pos.y = mark_pos[1]
            pos.z = mark_pos[2]
    if DEBUG:
        print(f'LOCAL WATER CHAINS: {moving_ligand_id} : {local_water_chains}')
        print(f'CHAIN NAME TO NAME: {moving_ligand_id} : {chain_name_to_chain}')

    # Drop residues and non-local waters from non-binding chains containing site waters
    new_chains = []
    for _model in new_structure:
        for _chain in _model:
            # Create a new chain to hold symmetry waters from non-binding chains
            new_chain = gemmi.Chain(_chain.name)
            # Iterate residues in the old chain, adding the local waters
            for _residue in _chain:
                if (_chain.name, str(_residue.seqid.num)+icode_to_string(_residue.seqid.icode)) in other_ligand_ids:
                    if DEBUG:
                        print(f'RESIDUE SKIP: {moving_ligand_id} : skipping {(_chain.name, str(_residue.seqid.num)+icode_to_string(_residue.seqid.icode))}')
                    continue
                if _residue.name == 'HOH':
                    if _chain.name in local_water_chains:
                        if str(_residue.seqid.num)+icode_to_string(_residue.seqid.icode) in local_water_chains[_chain.name]:
                            new_chain.add_residue(_residue.clone())
                else:
                    # Only add the chains that are part of the biological assemly the ligand
                    # is modelled as part of
                    if (_chain.name in lig_assembly.chains):
                        # If ligand drop other altlocs
                        if (_chain.name, str(_residue.seqid.num)+icode_to_string(_residue.seqid.icode)) == (moving_ligand_id[1], moving_ligand_id[2]):
                            new_residue = _residue.clone()
                            atoms_to_delete = set(
                                [
                                    (_atom.name, _atom.altloc) 
                                    for _atom in new_residue 
                                    if _atom.altloc != moving_ligand_id[3]
                                    ]
                            )
                            for _atom_name, _altloc in atoms_to_delete:
                                new_residue.remove_atom(_atom_name, _altloc)
                            # if len(atoms_to_delete) > 0:
                            #     print([len(_residue), len(new_residue), moving_ligand_id[3], atoms_to_delete])
                            #     exit()
                            new_chain.add_residue(new_residue)
                        else:
                            # Don't include other ligands
                            new_chain.add_residue(_residue.clone())
                    

            if len(new_chain) > 0:
                new_chains.append(new_chain)

    for new_chain in new_chains:
        del new_structure[0][new_chain.name]
        new_structure[0].add_chain(new_chain)

    if DEBUG:
        chain_ress = {_chain.name: len([x for x in _chain]) for _chain in new_structure[0]}
        print(f'CHAIN RESIUDES: {moving_ligand_id} : {chain_ress}')


    return new_structure


def get_alignment_transform(
    g, 
    moving_ligand_id, 
    reference_ligand_id,
    neighbourhood_transforms,
    chain_to_assembly_transform,
    assembly_transform,
):
    shortest_path: list[tuple[str, str, str]] = nx.shortest_path(g, moving_ligand_id, reference_ligand_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = moving_ligand_id
    running_transform = gemmi.Transform()
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = transform_to_gemmi(
                neighbourhood_transforms[
                    (
                        next_ligand_id,
                        previous_ligand_id,
                    )
                ]
            )
            running_transform = transform.combine(running_transform)

        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Update the transform with the assembly alignment
    # # Get the xtalform to assembly transform
    running_transform = alignment_heirarchy._transform_to_gemmi(chain_to_assembly_transform).combine(running_transform)

    # # Get the assembly alignment transform
    running_transform = alignment_heirarchy._transform_to_gemmi(assembly_transform).combine(running_transform)

    return running_transform


def _align_structure(
    _structure,
    moving_ligand_id: tuple[str, str, str],
    reference_ligand_id: tuple[str, str, str],
    neighbourhood: dt.Neighbourhood,
    dataset_ligand_neighbourhood_ids,
    g,
    neighbourhood_transforms: dict[tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform],
    # conformer_site_transforms: dict[tuple[str, str], dt.Transform],
    # canonical_site_id: str,
    # conformer_site_id: str,
    xtalform: dt.XtalForm,
    out_path: Path,
    # site_reference_xform,
    chain_to_assembly_transform,
    assembly_transform,
    xtalform_sites,
):


    # Drop chains that are not in the ligand's biological assembly
    reduced_structure = _drop_non_assembly_chains_and_symmetrize_waters(
        _structure, neighbourhood, moving_ligand_id, dataset_ligand_neighbourhood_ids, xtalform, xtalform_sites
    )

    transform = get_alignment_transform(
        g, 
        moving_ligand_id, 
        reference_ligand_id,
        neighbourhood_transforms,
        chain_to_assembly_transform,
        assembly_transform,
    )

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(transform)}")

    aligned_structure = superpose_structure(transform, reduced_structure)

    # Write the fully aligned structure
    aligned_structure.write_pdb(str(out_path))


def _align_reference_structure(
    _structure,
    dtag: str,
    reference_structure_transforms: dict[tuple[str, str], dt.Transform],
    canonical_site_id: str,
    out_path: Path,
):
    running_transform = transform_to_gemmi(reference_structure_transforms[(dtag, canonical_site_id)])

    # Site alignment transform
    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(running_transform)}")

    new_structure = superpose_structure(running_transform, _structure)

    # Write the fully aligned structure
    new_structure.write_pdb(str(out_path))


def _align_structures_from_sites(
    structures,
    canonical_sites: CanonicalSites,
    conformer_sites: ConformerSites,
    transforms: Transforms,
    neighbourhoods: LigandNeighbourhoods,
    xtalforms: XtalForms,
    assigned_xtalforms: AssignedXtalForms,
    g,
    site_transforms: SiteTransforms,
    output: Output,
):
    # Iterate sites
    for canonical_site_id, canonical_site in canonical_sites.iter():
        logger.debug(f"Canonical Site id is: {canonical_site_id}")

        for conformer_site_id, conformer_site in conformer_sites.iter():
            if conformer_site_id not in canonical_site.subsite_ids:
                continue
            logger.debug(f"Conformer Site id is: {canonical_site_id}")

            ligand_ids = conformer_site.members
            logger.debug(f"Conformer Site members: {len(conformer_site.members)}")

            # Select the alignment reference ligand_id
            # TODO: Track reference properly
            reference_ligand_id = conformer_site.reference_ligand_id

            # For each other ligand
            for moving_ligand_id in ligand_ids:
                logger.info(f"Alligning ligand: {moving_ligand_id}")
                # Get the shortest alignment path to the reference
                # Initial structure
                _structure = structures[moving_ligand_id.dtag].clone()

                # Get output path
                aod = Path(output.source_dir)
                output_path = (
                    aod
                    / output.dataset_output[moving_ligand_id.dtag][moving_ligand_id.chain][
                        moving_ligand_id.residue
                    ].aligned_structures[canonical_site_id]
                )
                # Align the ligand
                align_structure(
                    _structure,
                    moving_ligand_id,
                    reference_ligand_id,
                    g,
                    transforms,
                    site_transforms,
                    canonical_site_id,
                    conformer_site_id,
                    output_path,
                )


def get_structure_from_template(structure):
    new_structure = structure.clone()
    chain_names = [chain.name for chain in new_structure[0]]
    model = new_structure[0]
    for chain_name in chain_names:
        del model[chain_name]
    
    num_chains = len([x for x in new_structure[0]])
    assert num_chains == 0, f'Should have no chains in template structure! Got {num_chains}'

    return new_structure


def transform_chain(structure, chain_name, image):
    new_chain = gemmi.Chain(chain_name)

    transform = dt.transform_to_gemmi(image)

    for residue in structure[0][chain_name]:
        new_residue = residue.clone()
        for atom in new_residue:
            pos = atom.pos
            pos_frac = structure.cell.fractionalize(pos)
            pos_frac_transformed = gemmi.Fractional(transform.apply(pos_frac))
            new_pos = structure.cell.orthogonalize(pos_frac_transformed)
            atom.pos = new_pos
        new_chain.add_residue(new_residue)

    return new_chain


def get_structure_from_chain_images(structure, artefact_chains):
    # Get a clean template structure
    new_structure = get_structure_from_template(structure)
    
    # Build the transformed chains
    new_chains = []
    chain_counts = {chain_name: 0 for chain_name, _ in artefact_chains}
    for (chain_name, image_name), image in artefact_chains.items():
        transformed_chain = transform_chain(structure, chain_name, image)
        transformed_chain.name = f'{chain_name}{chain_counts[chain_name]}'
        new_chains.append(transformed_chain)
        
        chain_counts[chain_name] += 1

    # Add the transformed chains with unique names 
    for chain in new_chains:
        new_structure[0].add_chain(chain)

    num_chains = len([x for x in new_structure[0]])
    assert num_chains == len(artefact_chains), f'Should have {len(artefact_chains)} but have {num_chains}'

    return new_structure



def _align_artefacts(
    _structure,
    moving_ligand_id: tuple[str, str, str],
    reference_ligand_id: tuple[str, str, str],
    neighbourhood: dt.Neighbourhood,
    dataset_ligand_neighbourhood_ids,
    g,
    neighbourhood_transforms: dict[tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform],
    xtalform: dt.XtalForm,
    out_path: Path,
    # site_reference_xform,
    chain_to_assembly_transform,
    assembly_transform,
    xtalform_sites,
):
    # Get the artefact chains in the neighbourhood
    identity_transform = dt.Transform(
        [0.0,0.0,0.0],
        [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    )
    identity_transform_string = dt.transform_to_string(identity_transform)
    artefact_chains = {
        (atom_id[0], dt.transform_to_string(atom.image)): atom.image 
        for atom_id, atom 
        in neighbourhood.artefact_atoms.items()
        if dt.transform_to_string(atom.image) != identity_transform_string
        }
    # print(f'{moving_ligand_id[0]} : {artefact_chains}')

    # Build the artefact structure
    artefact_structure = get_structure_from_chain_images(
        _structure, 
        artefact_chains,
        )

    # Align it with the same transform as was used for the non-artefact atoms
    transform = get_alignment_transform(
        g, 
        moving_ligand_id, 
        reference_ligand_id,
        neighbourhood_transforms,
        chain_to_assembly_transform,
        assembly_transform,
    )

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(transform)}")

    aligned_structure = superpose_structure(
        transform, 
        artefact_structure,
        )

    # Write the fully aligned structure
    aligned_structure.write_pdb(str(out_path))
