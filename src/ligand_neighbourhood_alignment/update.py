import sys
from pathlib import Path

from loguru import logger
from rich import print
from joblib import Parallel, delayed

# in the original code, the imports were like this:
# from rich import print as rprint

# logger.remove()  # for someone not familiar with the lib, whats going on here?
# logger.add(sys.stdout, level="WARNING")
# from rich import print
# I wonder if this has some significance?

import gemmi

from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.forced_alignments import (
    _update_reference_structure_transforms,
)
from ligand_neighbourhood_alignment import alignment_landmarks
from ligand_neighbourhood_alignment.structure_alignment import (
    _align_structure,
    _align_reference_structure,
    _align_artefacts,
)
from ligand_neighbourhood_alignment.map_alignment import read_xmap, read_xmap_from_mtz, __align_xmap
from ligand_neighbourhood_alignment import alignment_heirarchy
from ligand_neighbourhood_alignment.io import (
    _get_structures,
    _save_assignments,
    _save_neighbourhoods,
    _save_ligand_neighbourhood_transforms,
    _save_graph,
    _save_connected_components,
    _save_conformer_sites,
    _save_canonical_sites,
    _save_xtalform_sites,
    _save_reference_structure_transforms,
    _save_fs_model,
    save_yaml,
)
from ligand_neighbourhood_alignment.ligand_neighbourhoods import (
    _get_neighbourhoods,
    _update_ligand_neighbourhood_transforms,
)
from ligand_neighbourhood_alignment.conformer_sites import _update_conformer_sites, _get_connected_components
from ligand_neighbourhood_alignment.xtalform_sites import _update_xtalform_sites
from ligand_neighbourhood_alignment.canonical_sites import _update_canonical_sites
from ligand_neighbourhood_alignment.xtalform_assignment import _assign_dataset
from ligand_neighbourhood_alignment.neighbourhood_graph import _update_graph
from ligand_neighbourhood_alignment.fs import _update_fs_model

logger.remove()  # for someone not familiar with the lib, whats going on here?
logger.add(sys.stdout, level="WARNING")


def perform_all_alignments(
    dtag,
    chain,
    residue,
    altloc,
    version,
    ligand_neighbourhood_output,
    canonical_site_id,
    fs_model,
    structures,
    conformer_sites,
    canonical_sites,
    dataset_assignments,
    xtalforms,
    xtalform_sites,
    ligand_neighbourhoods,
    alignability_graph,
    ligand_neighbourhood_transforms,
    chain_to_assembly_transforms,
    assembly_transforms,
    datasets,
    reference_xmap_path,
):
    reference_xmap = read_xmap_from_mtz(reference_xmap_path)

    aligned_structure_path = ligand_neighbourhood_output.aligned_structures[canonical_site_id]
    aligned_artefacts_path = ligand_neighbourhood_output.aligned_artefacts[canonical_site_id]
    aligned_event_map_path = ligand_neighbourhood_output.aligned_event_maps[canonical_site_id]

    if not (fs_model.source_dir.parent / aligned_structure_path).exists() | Path(aligned_structure_path).exists():
        _structure = structures[dtag].clone()
        canonical_site = canonical_sites[canonical_site_id]

        # Check for the matching conformer site
        conformer_site = None
        for conformer_site_id in canonical_site.conformer_site_ids:
            if (dtag, chain, residue, altloc, version) in conformer_sites[conformer_site_id].members:
                conformer_site = conformer_sites[conformer_site_id]
                break
        if conformer_site is None:
            print(f"Skipping alignment of {dtag} {chain} {residue} {altloc} to site {canonical_site_id}!")
            return
        moving_ligand_id = (dtag, chain, residue, altloc, version)
        reference_ligand_id = conformer_site.reference_ligand_id

        # Get the site chain
        site_chain = None
        xtalform_site = None
        site_reference_ligand_id = conformer_site.reference_ligand_id
        site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
        site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]
        for xsid, _xtalform_site in xtalform_sites.items():
            _xtalform_id = _xtalform_site.xtalform_id
            _xtalform_canonical_site_id = _xtalform_site.canonical_site_id
            if (
                (_xtalform_id == site_reference_ligand_xtalform_id)
                & (_xtalform_canonical_site_id == canonical_site_id)
                & (site_reference_ligand_id in _xtalform_site.members)
            ):
                xtalform_site = _xtalform_site
        site_chain = xtalform_site.crystallographic_chain

        # Aligns to conformer site, then to the corresponding assembly, then from that assembly to
        # the global frame
        try:
            _align_structure(
                _structure,
                moving_ligand_id,
                reference_ligand_id,
                ligand_neighbourhoods[moving_ligand_id],
                [nid for nid in ligand_neighbourhoods if nid[0] == dtag],
                alignability_graph,
                ligand_neighbourhood_transforms,
                xtalforms[dataset_assignments[dtag]],
                aligned_structure_path,
                chain_to_assembly_transform=chain_to_assembly_transforms[
                    (conformer_site.reference_ligand_id[0], site_chain)
                ],
                assembly_transform=assembly_transforms[
                    xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                    .assemblies[
                        alignment_heirarchy._chain_to_xtalform_assembly(
                            site_chain,
                            xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                        )
                    ]
                    .assembly
                ],
                xtalform_sites=xtalform_sites,
            )

        except:
            logger.info(f"Failed to generate aligned structure {aligned_structure_path}")
            raise

        # Perform the map alignments

        _align_artefacts(
            _structure,
            moving_ligand_id,
            reference_ligand_id,
            ligand_neighbourhoods[moving_ligand_id],
            [nid for nid in ligand_neighbourhoods if nid[0] == dtag],
            alignability_graph,
            ligand_neighbourhood_transforms,
            xtalforms[dataset_assignments[dtag]],
            aligned_artefacts_path,
            chain_to_assembly_transform=chain_to_assembly_transforms[
                (conformer_site.reference_ligand_id[0], site_chain)
            ],
            assembly_transform=assembly_transforms[
                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                .assemblies[
                    alignment_heirarchy._chain_to_xtalform_assembly(
                        site_chain,
                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                    )
                ]
                .assembly
            ],
            xtalform_sites=xtalform_sites,
        )
        # else:
        #     logger.info(f"Already output aligned artefacts!")
        # else:
        #     logger.info(f"Already output structure!")

        # Align the xmaps
        try:
            st_path = fs_model.source_dir.parent / aligned_structure_path
            if not st_path.exists():
                st_path = aligned_structure_path
            aligned_structure = gemmi.read_structure(str(st_path))
            aligned_res = aligned_structure[0][chain][str(residue)][0]
            xmap_path = datasets[dtag].ligand_binding_events[(dtag, chain, residue, altloc)].xmap
        except Exception as e:
            print(f'Structure path is: {st_path}')
            print(f'canonical site id is: {canonical_site_id}')
            print(f'ligand binding events are: {[x for x in datasets[dtag].ligand_binding_events.keys()]}')
            print(f'Error occured for: {dtag} {chain} {residue} {altloc}')
            print(f'Structure chains are: {[_chain.name for _chain in aligned_structure[0]]}')
            print(f'Chain residues are: {[_res.seqid.num for _res in aligned_structure[0][chain]]}')
            # print(f'Cha[])
            raise e

        if (xmap_path != "None") and (xmap_path is not None):
            xmap = read_xmap(xmap_path)

            __align_xmap(
                ligand_neighbourhoods[(dtag, chain, residue, altloc, version)],
                alignability_graph,
                ligand_neighbourhood_transforms,
                reference_xmap,
                reference_ligand_id,
                moving_ligand_id,
                xmap,
                aligned_event_map_path,
                ligand_neighbourhood_output.aligned_event_maps_crystallographic[canonical_site_id],
                aligned_res,
                chain_to_assembly_transform=chain_to_assembly_transforms[
                    (conformer_site.reference_ligand_id[0], site_chain)
                ],
                assembly_transform=assembly_transforms[
                    xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                    .assemblies[
                        alignment_heirarchy._chain_to_xtalform_assembly(
                            site_chain,
                            xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                        )
                    ]
                    .assembly
                ],
            )
            mtz_path = datasets[dtag].mtz
            if mtz_path and mtz_path != "None":
                xmap = read_xmap_from_mtz(mtz_path, "2Fo-Fc")
                __align_xmap(
                    ligand_neighbourhoods[(dtag, chain, residue, altloc, version)],
                    alignability_graph,
                    ligand_neighbourhood_transforms,
                    reference_xmap,
                    reference_ligand_id,
                    moving_ligand_id,
                    xmap,
                    ligand_neighbourhood_output.aligned_xmaps[canonical_site_id],
                    ligand_neighbourhood_output.aligned_xmaps_crystallographic[canonical_site_id],
                    aligned_res,
                    chain_to_assembly_transform=chain_to_assembly_transforms[
                        (conformer_site.reference_ligand_id[0], site_chain)
                    ],
                    assembly_transform=assembly_transforms[
                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                        .assemblies[
                            alignment_heirarchy._chain_to_xtalform_assembly(
                                site_chain,
                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                            )
                        ]
                        .assembly
                    ],
                )
                xmap = read_xmap_from_mtz(mtz_path, "Fo-Fc")
                __align_xmap(
                    ligand_neighbourhoods[(dtag, chain, residue, altloc, version)],
                    alignability_graph,
                    ligand_neighbourhood_transforms,
                    reference_xmap,
                    reference_ligand_id,
                    moving_ligand_id,
                    xmap,
                    ligand_neighbourhood_output.aligned_diff_maps[canonical_site_id],
                    ligand_neighbourhood_output.aligned_diff_maps_crystallographic[canonical_site_id],
                    aligned_res,
                    chain_to_assembly_transform=chain_to_assembly_transforms[
                        (conformer_site.reference_ligand_id[0], site_chain)
                    ],
                    assembly_transform=assembly_transforms[
                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                        .assemblies[
                            alignment_heirarchy._chain_to_xtalform_assembly(
                                site_chain,
                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                            )
                        ]
                        .assembly
                    ],
                )

        else:
            logger.info(f"Already output xmap!")


def update(
    fs_model: dt.FSModel,
    datasets: dict[str, dt.Dataset],
    reference_datasets: dict[str, dt.Dataset],
    new_datasets: dict[str, dt.Dataset],
    assemblies: dict[str, dt.Assembly],
    xtalforms: dict[str, dt.XtalForm],
    dataset_assignments: dict[str, str],
    ligand_neighbourhoods: dict[dt.LigandNeighbourhoodID, dt.Neighbourhood],
    alignability_graph,
    connected_components,
    ligand_neighbourhood_transforms: dict[tuple[dt.LigandNeighbourhoodID, dt.LigandNeighbourhoodID], dt.Transform],
    conformer_sites: dict[str, dt.ConformerSite],
    canonical_sites: dict[str, dt.CanonicalSite],
    xtalform_sites: dict[str, dt.XtalFormSite],
    reference_structure_transforms: dict[tuple[str, str], dt.Transform],
    assembly_landmarks,
    assembly_transforms,
    version,
    debug=False,
):
    """
    The main interface for XCA to LNA functionality.

    This function performs the crystallographic functions of an update step in XCA, including:
    1.  Deriving the alignment hierarchy. This creates an artifical "reference assembly" that defines
        the relative positions of biological chains to one another in the final alignment.
        Note that the binding site alignments will then be relative to this articial assembly
        of chains, maximizing the -local- goodness of fit to the binding sites
    2.  Generate the landmarks that will be used for alignment.
    3.  Determine how to transform each assembly onto the reference assembly. This way binding sites
        can be aligned to their own biological assembly and then to the reference assembly.
    4.  Assign crystalforms to datasets. This is a heuristic to save people manually assigning
        a crystalform. It seems fairly reliable.
    5.  Get the neighbourhoods around the ligand. This will determine both the local biological
        residues and local crystalographic residues.
    6.  Get the chain to assembly transforms. These will map each chain in each dataset to their
        local assembly (which can then be mapped to the reference assembly).
    7.  Determine the transforms between ligand neighbourhoods and by extension the graph of
        which transforms are possible.
    8.  Get the conformer sites
    9.  Get the canonical sites
    10. Get the xtalform sites
    11. Get the transforms from forced ligand-less datasets
    12. Update the fs model
    13. Generate the alignments of structures and xmaps.
    """

    logger.info(f"Version is: {version}")
    # Get the structures
    structures: dict = _get_structures(datasets)

    # Get the assembly alignment hierarchy
    hierarchy, biochain_priorities = alignment_heirarchy._derive_alignment_heirarchy(assemblies)
    save_yaml(fs_model.hierarchy, hierarchy, lambda x: x)
    save_yaml(fs_model.biochain_priorities, biochain_priorities, lambda x: x)

    # Get the assembly hierarchy transforms
    for assembly_name, assembly in assemblies.items():
        # Do not update if already have landmarks!
        if assembly_name in assembly_landmarks:
            continue
        as_st = alignment_heirarchy._get_assembly_st(assembly, structures[assembly.reference])
        assembly_landmarks[assembly_name] = alignment_landmarks.structure_to_landmarks(as_st)

    if debug:
        print(assembly_landmarks)

    save_yaml(fs_model.assembly_landmarks, assembly_landmarks, alignment_landmarks.assembly_landmarks_to_dict)

    for assembly_name, assembly in assemblies.items():
        # Do not update if already have assembly transform!
        if assembly_name in assembly_transforms:
            continue
        assembly_transforms[assembly_name] = alignment_heirarchy._calculate_assembly_transform_sequence(
            hierarchy,
            assembly_name,
            assembly_landmarks,
        )
    save_yaml(fs_model.assembly_transforms, assembly_transforms, lambda x: x)

    # Assign datasets
    new_dataset_assignments = {}
    for dtag, dataset in new_datasets.items():
        new_dataset_assignments[dtag] = _assign_dataset(
            dataset,
            assemblies,
            xtalforms,
            structures[dtag],
            structures,
        )
    dataset_assignments.update(new_dataset_assignments)

    _save_assignments(fs_model, dataset_assignments)
    logger.info(f"Assigned {len(dataset_assignments)} xtalform assignments to datasets!")

    # Get neighbourhoods
    logger.info(f"Updating neighbourhoods")
    for dtag, dataset in new_datasets.items():
        if debug:
            print(f'\tGetting neighbourhood: {dtag}')
        xtalform = xtalforms[dataset_assignments[dtag]]
        neighborhoods = _get_neighbourhoods(dataset, xtalform, assemblies, version)
        logger.info(f"Dataset {dtag} has {len(neighborhoods)} ligand neighbourhoods")
        for lid, neighbourhood in neighborhoods.items():
            ligand_neighbourhoods[lid] = neighbourhood
    logger.info(f"Found {len(ligand_neighbourhoods)} ligand neighbourhoods!")
    _save_neighbourhoods(fs_model, ligand_neighbourhoods)

    if debug:
        print('Assembly Landmarks')
        print(assembly_landmarks)

    # Get chain to assembly transforms
    logger.info(f"Getting chain-to-assembly transforms...")
    chain_to_assembly_transforms = {}
    for dtag, st in structures.items():
        xtalform_chains = [
            _chain
            for _xassembly in xtalforms[dataset_assignments[dtag]].assemblies.values()
            for _chain in _xassembly.chains
        ]

        dataset_chains = [_chain.name for _chain in st[0]]
        dataset_ligand_chains = [_x[1] for _x in ligand_neighbourhoods if _x[0] == dtag]
        for _chain in xtalform_chains:
            if _chain not in xtalform_chains:
                raise Exception(
                    f"A xtalform assignment error has occured. Dataset {dtag} has chain {_chain} in its chains {dataset_chains} however its assigned xtalform {dataset_assignments[dtag]} has chain {xtalform_chains}"
                )
            try:
                if debug & (dtag == 'A71EV2A-x7510'):
                    _debug = True
                else:
                    _debug = False

                print('dtag')
                print(dtag)

                chain_to_assembly_transforms[
                    (
                        dtag,
                        _chain,
                    )
                ] = alignment_heirarchy._get_structure_chain_to_assembly_transform(
                    st, _chain, xtalforms[dataset_assignments[dtag]], assemblies, assembly_landmarks, _debug
                )
            except Exception as e:
                print(f'Exception in dataset: {dtag}, in xtalform {dataset_assignments[dtag]}')
                raise e
    logger.info(f'Got {len(chain_to_assembly_transforms)} chain to assembly transforms')
    save_yaml(
        fs_model.chain_to_assembly,
        chain_to_assembly_transforms,
        alignment_heirarchy.chain_to_assembly_transforms_to_dict,
    )
    logger.info(f'Got {len(chain_to_assembly_transforms)} chain to assembly transforms')

    # Update graph
    logger.info(f"Updating alignment graph...")
    logger.info(f"Previously had {len(ligand_neighbourhood_transforms)} alignments between neighbourhoods")
    for lid, neighbourhood in ligand_neighbourhoods.items():
        _update_ligand_neighbourhood_transforms(
            ligand_neighbourhood_transforms,
            lid,
            ligand_neighbourhoods,
            structures,
        )

    logger.info(f"Now have {len(ligand_neighbourhood_transforms)} alignments between neighbourhoods")
    _save_ligand_neighbourhood_transforms(fs_model, ligand_neighbourhood_transforms)

    # Update the alignment graph
    logger.info(f"Updating alignment graph...")
    logger.info(f"Previously had {len(alignability_graph.nodes)} nodes")
    logger.info(f"Previously had {len(alignability_graph.edges)} edges")
    _update_graph(
        alignability_graph,
        ligand_neighbourhoods,
        ligand_neighbourhood_transforms,
    )
    logger.info(f"Now have {len(alignability_graph.nodes)} nodes")
    logger.info(f"Now have {len(alignability_graph.edges)} edges")
    _save_graph(fs_model, alignability_graph)

    # Update conformer sites
    logger.info(f"Updating conformer sites...")
    connected_components = _get_connected_components(alignability_graph, connected_components)
    _save_connected_components(fs_model, connected_components)
    logger.info(f"Got {len(connected_components)} connected components")
    logger.info(f"Previously had {len(conformer_sites)} conformer sites")

    for connected_component_id, connected_component in connected_components.items():
        # Update new datasets to indicate everything sharing a connected component
        # Match new component to old ones by membership, and expand old ones if available otherwise create new one
        _update_conformer_sites(
            conformer_sites,
            connected_component_id,
            connected_component,
            ligand_neighbourhoods,
            structures,
            assemblies,
            xtalforms,
            dataset_assignments,
        )
    logger.info(f"Now have {len(conformer_sites)} conformer sites")
    _save_conformer_sites(fs_model, conformer_sites)

    # Update canonical sites
    logger.info(f"Previously had {len(canonical_sites)} canonical sites")
    for conformer_site_id, conformer_site in conformer_sites.items():
        # If conformer site in a canonical site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        _update_canonical_sites(canonical_sites, conformer_site, conformer_site_id, ligand_neighbourhoods)
    logger.info(f"Now have {len(canonical_sites)} canonical sites")
    logger.info(f"Global reference dtag is: {list(canonical_sites.values())[0].global_reference_dtag}")
    _save_canonical_sites(fs_model, canonical_sites)

    # Update crystalform sites
    logger.info(f"Previously had {len(xtalform_sites)} xtalform sites")
    for canonical_site_id, canonical_site in canonical_sites.items():
        # If canonical site in a xtalform site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        _update_xtalform_sites(
            xtalform_sites,
            canonical_site,
            canonical_site_id,
            dataset_assignments,
            conformer_sites,
            ligand_neighbourhoods,
            debug,
        )
    logger.info(f"Now have {len(xtalform_sites)} xtalform sites")
    _save_xtalform_sites(fs_model, xtalform_sites)

    # Update the reference structure transforms
    logger.info(f"Previously had {len(reference_structure_transforms)} reference structure transforms")
    for dtag, dataset in reference_datasets.items():
        for canonical_site_id, canonical_site in canonical_sites.items():
            key = (dtag, canonical_site_id)
            if key not in reference_structure_transforms:
                _update_reference_structure_transforms(
                    reference_structure_transforms,
                    key,
                    structures,
                    canonical_site,
                    conformer_sites,
                    assemblies,
                    xtalforms,
                    dataset_assignments,
                    xtalform_sites,
                    canonical_site_id,
                )
    logger.info(f"Now have {len(reference_structure_transforms)} reference structure transforms")
    _save_reference_structure_transforms(
        fs_model,
        reference_structure_transforms,
    )

    # Update output: check if aligned data for each lid in canon site is already there and if not add it
    _update_fs_model(fs_model, canonical_sites, conformer_sites, reference_datasets, reference_structure_transforms)
    _save_fs_model(fs_model)

    # Perform the alignments
    reference_xmap_path = datasets[[x for x in canonical_sites.values()][0].global_reference_dtag].mtz

    fs = []
    for dtag, dataset_alignment_info in fs_model.alignments.items():
        for chain, chain_alignment_info in dataset_alignment_info.items():
            for residue, residue_alignment_info in chain_alignment_info.items():
                for altloc, altloc_alignment_info in residue_alignment_info.items():
                    for version, ligand_neighbourhood_output in altloc_alignment_info.items():
                        for (
                            canonical_site_id,
                            aligned_structure_path,
                        ) in ligand_neighbourhood_output.aligned_structures.items():
                            fs.append(
                                delayed(perform_all_alignments)(
                                    dtag,
                                    chain,
                                    residue,
                                    altloc,
                                    version,
                                    ligand_neighbourhood_output,
                                    canonical_site_id,
                                    fs_model,
                                    structures,
                                    conformer_sites,
                                    canonical_sites,
                                    dataset_assignments,
                                    xtalforms,
                                    xtalform_sites,
                                    ligand_neighbourhoods,
                                    alignability_graph,
                                    ligand_neighbourhood_transforms,
                                    chain_to_assembly_transforms,
                                    assembly_transforms,
                                    datasets,
                                    reference_xmap_path,
                                )
                            )

    Parallel(verbose=100, n_jobs=-1)(f for f in fs)

    # Generate alignments of references to each canonical site
    for dtag, dataset_alignment_info in fs_model.reference_alignments.items():
        for canonical_site_id, alignment_info in dataset_alignment_info.items():
            aligned_structure_path = alignment_info["aligned_structures"]
            logger.info(f"Outputting reference structure: {aligned_structure_path}")
            if not (
                (fs_model.source_dir.parent / aligned_structure_path).exists() | Path(aligned_structure_path).exists()
            ):
                _structure = structures[dtag].clone()
                _align_reference_structure(
                    _structure,
                    dtag,
                    reference_structure_transforms,
                    canonical_site_id,
                    alignment_info["aligned_structures"],
                )
            else:
                logger.info(f"Already output reference structure!")

    # # Align the maps

    # logger.info(f"Outputting xmaps...")
    # for dtag, dataset_alignment_info in fs_model.alignments.items():
    #     for chain, chain_alignment_info in dataset_alignment_info.items():
    #         for residue, residue_alignment_info in chain_alignment_info.items():
    #             for version, ligand_neighbourhood_output in residue_alignment_info.items():
    #                 for (
    #                     canonical_site_id,
    #                     aligned_event_map_path,
    #                 ) in ligand_neighbourhood_output.aligned_event_maps.items():
    #                     logger.info(f"Writing to: {aligned_event_map_path}")

    #                     rel = True
    #                     try:
    #                         rel = Path(aligned_event_map_path).is_relative_to(fs_model.source_dir)
    #                     except:
    #                         rel = False
    #                     if rel:
    #                         _structure = structures[dtag].clone()
    #                         canonical_site = canonical_sites[canonical_site_id]
    #                         # Check for the matching conformer site
    #                         conformer_site = None
    #                         for conformer_site_id in canonical_site.conformer_site_ids:
    #                             if (dtag, chain, residue, version) in conformer_sites[conformer_site_id].members:
    #                                 conformer_site = conformer_sites[conformer_site_id]
    #                                 break

    #                         if conformer_site is None:
    #                             print(f"Skipping alignment of {dtag} {chain} {residue} to site {canonical_site_id}!")
    #                             continue

    #                         moving_ligand_id = (dtag, chain, residue, version)
    #                         reference_ligand_id = conformer_site.reference_ligand_id

    #                         xmap_path = datasets[dtag].ligand_binding_events[(dtag, chain, residue)].xmap

    #                         aligned_structure_path = ligand_neighbourhood_output.aligned_structures[canonical_site_id]
    #                         st_path = fs_model.source_dir.parent / aligned_structure_path
    #                         if not st_path.exists():
    #                             st_path = aligned_structure_path
    #                         aligned_structure = gemmi.read_structure(str(st_path))
    #                         aligned_res = aligned_structure[0][chain][str(residue)][0]

    #                         # Get the site chain
    #                         site_chain = None
    #                         xtalform_site = None
    #                         site_reference_ligand_id = conformer_site.reference_ligand_id
    #                         site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
    #                         site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]
    #                         for xsid, _xtalform_site in xtalform_sites.items():
    #                             _xtalform_id = _xtalform_site.xtalform_id
    #                             _xtalform_canonical_site_id = _xtalform_site.canonical_site_id
    #                             if (
    #                                 (_xtalform_id == site_reference_ligand_xtalform_id)
    #                                 & (_xtalform_canonical_site_id == canonical_site_id)
    #                                 & (site_reference_ligand_id in _xtalform_site.members)
    #                             ):
    #                                 xtalform_site = _xtalform_site
    #                         site_chain = xtalform_site.crystallographic_chain

    return fs_model
