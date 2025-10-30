from loguru import logger

from ligand_neighbourhood_alignment.data import (
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SystemSites,
)
from src.ligand_neighbourhood_alignment.alignment_core import (
    match_neighbourhood_to_site,
    match_neighbourhood_to_sites,
)


def get_canonical_sites(
    initial_system_sites: SystemSites | None,
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood],
) -> dict[int, CanonicalSite]:
    # Get the existing canonical sites
    canonical_sites: dict[int, CanonicalSite] = {}
    if initial_system_sites:
        canonical_sites = initial_system_sites.canonical_site

    # Get the number of canonical sites
    if len(canonical_sites) == 0:
        canonical_site_num = 0
    else:
        canonical_site_num = max(canonical_sites)
    logger.info(f"Number of canon sites: {canonical_site_num}")

    # Iterate neighbourhoods, matching them to existing sites, and
    # creating a new site if necessary
    canonical_site_members: dict[int, list[LigandID]] = {}
    for ligand_id, ligand_neighbourhood in ligand_neighbourhoods.items():
        logger.debug(f"{ligand_id}")
        # Check if there is a match
        match: int | None = match_neighbourhood_to_sites(canonical_sites, ligand_neighbourhood)

        # If so add the ligand id to the members of the canonical site
        if match:
            # Check that there are already matches, and if not make the
            # list to add them to
            if match not in canonical_site_members:
                canonical_site_members[match] = []
            canonical_site_members[match].append(ligand_id)
            canonical_sites[match].members.append(ligand_id)

        # Otherwise create a new canonical site templated on the ligand
        # neighbourhood
        else:
            canonical_sites[canonical_site_num] = CanonicalSite(
                id=canonical_site_num,
                name="",
                refpdb="",
                atoms=ligand_neighbourhood.atoms,
                literatureref="",
                members=[
                    ligand_id,
                ],
            )
            canonical_site_members[canonical_site_num] = [
                ligand_id,
            ]
            canonical_site_num += 1

    # Reiterate neighbourhoods, adding them to any extant sites that they
    # can match
    for ligand_id, ligand_neighbourhood in ligand_neighbourhoods.items():
        for (
            canonical_site_index,
            canonical_site_ligand_ids,
        ) in canonical_site_members.items():
            canonical_site: CanonicalSite = canonical_sites[canonical_site_index]
            rematch: bool = match_neighbourhood_to_site(canonical_site, ligand_neighbourhood)
            if rematch:
                if ligand_id not in canonical_site_ligand_ids:
                    canonical_site_ligand_ids.append(ligand_id)
                    canonical_site.members.append(ligand_id)

    return canonical_sites
