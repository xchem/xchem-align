# from loguru import logger

from ligand_neighbourhood_alignment.data import (  # XtalFormSite,
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SystemSites,
)

# from ligand_neighbourhood_alignment.matching import match_atoms


def get_xtal_form_sites(
    initial_system_sites: SystemSites | None,
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood],
    canonical_sites: dict[int, CanonicalSite],
):  # -> dict[int, XtalFormSite]:
    ...


#     # Get the current xtalform sites
#     xtal_form_sites: dict[int, XtalFormSite] = {}
#     if initial_system_sites:
#         xtal_form_sites = initial_system_sites.xtal_form_site

#     # Get the number of xtalform sites
#     if len(xtal_form_sites) == 0:
#         xtal_form_sites_num = 0
#     else:
#         xtal_form_sites_num = max(canonical_sites)
#     logger.info(f"Current number of xtalformsites: {xtal_form_sites_num}")

#     # Iterate through canonical sites, partitioning their members
# according to
#     # whether their artefact atoms can be aligned
#     for canonical_site_id, canonical_site in canonical_sites.items():
#         logger.debug(f"Processing canonical site: {canonical_site_id}")
#         # Determine whether each member aligns to an existing xtalform site
#         logger.debug(f"Num members: {len(canonical_site.members)}")
#         logger.debug(canonical_site.members)
#         for ligand_id in canonical_site.members:
#             ligand_neighbourhood = ligand_neighbourhoods[
#                 ligand_id
#             ]

#             # Check if a match to any xtalform site
#             match: bool = False
#             for xtal_form_id, xtal_form_site in xtal_form_sites.items():
#                 # TODO: Make this superposition of actual atoms
#                 match = match_atoms(
#                     xtal_form_site.artefact_atoms,
#                     ligand_neighbourhood.artefact_atoms,
#                 )
#                 if match:
#                     logger.debug(f"Site {ligand_id} matches {xtal_form_id}!")

#                     if ligand_id not in xtal_form_site.members:
#                         xtal_form_site.members.append(ligand_id)

#             # Otherwise create a new xtalform site and add the ligand to it
#             if not match:
#                 xtal_form_sites[xtal_form_sites_num] = XtalFormSite(
#                     id=xtal_form_sites_num,
#                     canon_site_id=canonical_site_id,
#                     xtal_form_id=0,
#                     code="",
#                     refpdb="",
#                     atoms=ligand_neighbourhood.atoms,
#                     artefact_atoms=ligand_neighbourhood.artefact_atoms,
#                     members=[
#                         ligand_id,
#                     ],
#                 )
#                 xtal_form_sites_num += 1

#     return xtal_form_sites
