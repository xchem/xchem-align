from ligand_neighbourhood_alignment.data import (
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SystemSites,
)


def get_xtal_form_sites(
    initial_system_sites: SystemSites | None,
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood],
    canonical_sites: dict[int, CanonicalSite],
):
    ...
