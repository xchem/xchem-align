from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment.ligand_neighbourhoods import _get_centroid_res


def _update_canonical_sites(
    canonical_sites: dict[str, dt.CanonicalSite],
    conformer_site: dt.ConformerSite,
    conformer_site_id,
    neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
):
    if len(canonical_sites) != 0:
        global_reference_dtag = [x for x in canonical_sites.values()][0].global_reference_dtag
    else:
        global_reference_dtag = conformer_site.reference_ligand_id[0]

    # If conformer site already in a canonical site skip
    for canonical_site_id, canonical_site in canonical_sites.items():
        if conformer_site_id in canonical_site.conformer_site_ids:
            return

            # Check each canonical site to see if conformer site already in it and if not
    # whether it shares enough residues to now be added
    matched = False
    conformer_site_residues = [(residue[1], residue[2]) for residue in conformer_site.residues]
    for canonical_site_id, canonical_site in canonical_sites.items():
        if matched:
            continue
        canonical_site_residues = [(residue[1], residue[2]) for residue in canonical_site.residues]
        if conformer_site_id not in canonical_site.conformer_site_ids:
            v = set(canonical_site_residues).intersection(set(conformer_site_residues))
            if len(v) >= 0.75 * len(canonical_site_residues):
                # Matched!
                matched = True
                canonical_site.conformer_site_ids.append(conformer_site_id)

    # If not matched to any existing canonical site create a new one
    if not matched:
        centroid_res = _get_centroid_res(conformer_site.residues, neighbourhoods[conformer_site.reference_ligand_id])
        canonical_site = dt.CanonicalSite(
            [
                conformer_site_id,
            ],
            conformer_site.residues,
            conformer_site_id,
            global_reference_dtag,
            (
                conformer_site.reference_ligand_id[0],
                centroid_res[0],
                centroid_res[1],
                conformer_site.reference_ligand_id[1],
            ),
        )

        canonical_site_id = conformer_site_id
        canonical_sites[canonical_site_id] = canonical_site
