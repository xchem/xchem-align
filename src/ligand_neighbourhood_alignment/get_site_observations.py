# import gemmi
# import loguru
# from ligand_neighbourhood_alignment.data import (
#     LigandID,
#     LigandNeighbourhood,
#     SiteObservation,
#     SystemSites,
#     XtalFormSite,
# )

# def get_site_observations(
#     system_sites: SystemSites | None,
#     xtal_form_sites: dict[int, XtalFormSite],
#     ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood],
# ) -> dict[LigandID, SiteObservation]:
#     # Get the current site observations

#     site_observations: dict[LigandID, SiteObservation] = {}
#     if system_sites:
#         site_observations = {lid: so for lid, so
#         in zip(system_sites.Sys ysystem_sites.site_observation

#     num_site_oservations: int = 0
#     if len(site_observations) != 0:
#         num_site_oservations = len(site_observations)

#     # For each member of each xtalformsite, create a new
#     for xtal_form_site_id, xtal_form_site in xtal_form_sites.items():
#         for ligand_id in xtal_form_site.members:
#             # If the site hasn't already been considered create a new
#             if ligand_id not in site_observations:
#                 # TODO: Fix dataset field
#                 site_observation: SiteObservation = SiteObservation(
#                     id=num_site_oservations,
#                     ligand_id=ligand_id,
#                     xtal_form_site_id=xtal_form_site_id,
#                     fragalysis_label="",
#                     description="",
#                     code="",
#                     dataset="",
#                 )
#                 site_observations[ligand_id] = site_observation
#                 num_site_oservations += 1

#     return site_observations
