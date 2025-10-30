from ligand_neighbourhood_alignment.data import LigandID


def merge_components(components: list[list[LigandID]]):
    for component in components:
        # Check against each new component and merge if overlaps heavily
        ...
