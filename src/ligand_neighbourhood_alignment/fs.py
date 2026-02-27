import os

from ligand_neighbourhood_alignment import dt, constants


def _update_fs_model(
    fs_model: dt.FSModel,
    canonical_sites: dict[str, dt.CanonicalSite],
    conformer_sites: dict[str, dt.ConformerSite],
    reference_datasets: dict[str, dt.Dataset],
    reference_transforms,
):
    # Iterate over canonical sites and their members, checking if they already have an output record and
    # if not creating one
    alignments = fs_model.alignments
    for canonical_site_id, canonical_site in canonical_sites.items():
        for conformer_site_id in canonical_site.conformer_site_ids:
            conformer_site = conformer_sites[conformer_site_id]
            for member in conformer_site.members:
                dtag, chain, residue, altloc, version = member
                if dtag not in alignments:
                    alignments[dtag] = {}
                if chain not in alignments[dtag]:
                    alignments[dtag][chain] = {}
                if residue not in alignments[dtag][chain]:
                    alignments[dtag][chain][residue] = {}
                if altloc not in alignments[dtag][chain][residue]:
                    alignments[dtag][chain][residue][altloc] = {}
                if version not in alignments[dtag][chain][residue][altloc]:
                    alignments[dtag][chain][residue][altloc][version] = dt.LigandNeighbourhoodOutput(
                        {}, {}, {}, {}, {}, {}, {}, {}
                    )

                ligand_neighbourhood_output: dt.LigandNeighbourhoodOutput = alignments[dtag][chain][residue][altloc][
                    version
                ]

                if not (fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                    os.mkdir(fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag)

                if canonical_site_id not in ligand_neighbourhood_output.aligned_structures:
                    ligand_neighbourhood_output.aligned_structures[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_artefacts[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_xmaps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_XMAP_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_diff_maps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_DIFF_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_event_maps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    # Crystallographic maps
                    ligand_neighbourhood_output.aligned_xmaps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_XMAP_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_diff_maps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_DIFF_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

                    ligand_neighbourhood_output.aligned_event_maps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_EVENT_MAP_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag,
                            chain=chain,
                            residue=residue,
                            altloc=dt.altloc_to_string(altloc),
                            version=version,
                            site=canonical_site_id,
                        )
                    )

    reference_alignments = fs_model.reference_alignments
    for dtag, dataset in reference_datasets.items():
        for canonical_site_id, canonical_site in canonical_sites.items():
            if (dtag, canonical_site_id) not in reference_transforms:
                continue

            if not (fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                os.mkdir(fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag)

            if dtag not in reference_alignments:
                reference_alignments[dtag] = {}

            if canonical_site_id not in reference_alignments[dtag]:
                reference_alignments[dtag][canonical_site_id] = {
                    "aligned_structures": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_STRUCTURE_TEMPLATE.format(dtag=dtag, site=canonical_site_id),
                    "aligned_artefacts": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                        dtag=dtag, site=canonical_site_id
                    ),
                    "aligned_xmaps": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_XMAP_TEMPLATE.format(dtag=dtag, site=canonical_site_id),
                }
