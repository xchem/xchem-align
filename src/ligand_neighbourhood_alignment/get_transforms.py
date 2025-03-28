import gemmi

from ligand_neighbourhood_alignment.data import (
    LigandID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Transform,
    Transforms,
)
from ligand_neighbourhood_alignment.matching import match_atom


def get_transform(
    reference_neighbourhood: LigandNeighbourhood,
    neighbourhood: LigandNeighbourhood,
):
    alignable_cas = {}
    for (
        ligand_1_atom_id,
        ligand_1_atom,
    ) in zip(reference_neighbourhood.atom_ids, reference_neighbourhood.atoms):
        for (
            ligand_2_atom_id,
            ligand_2_atom,
        ) in zip(neighbourhood.atom_ids, neighbourhood.atoms):
            if ligand_1_atom_id.atom == "CA":
                if match_atom(ligand_1_atom, ligand_2_atom, ignore_chain=True):
                    alignable_cas[ligand_1_atom_id] = (
                        gemmi.Position(
                            ligand_1_atom.x,
                            ligand_1_atom.y,
                            ligand_1_atom.z,
                        ),
                        gemmi.Position(
                            ligand_2_atom.x,
                            ligand_2_atom.y,
                            ligand_2_atom.z,
                        ),
                    )

    if len(alignable_cas) < 3:
        return gemmi.Transform(), []

    sup = gemmi.superpose_positions(
        [alignable_ca[0] for alignable_ca in alignable_cas.values()],
        [alignable_ca[1] for alignable_ca in alignable_cas.values()],
    )
    # logger.debug(f"Superposition: rmsd {sup.rmsd} n {len(alignable_cas)}")

    return sup.transform, [ligand_id for ligand_id in alignable_cas.keys()]


def get_transforms(ligand_neighbourhoods: LigandNeighbourhoods, g):
    transforms: dict[LigandID, dict[LigandID, Transform]] = {}
    for ligand_id_1, ligand_neighbourhood_1 in zip(
        ligand_neighbourhoods.ligand_ids,
        ligand_neighbourhoods.ligand_neighbourhoods,
    ):
        transforms[ligand_id_1] = {}
        for (
            ligand_id_2,
            ligand_neighbourhood_2,
        ) in zip(
            ligand_neighbourhoods.ligand_ids,
            ligand_neighbourhoods.ligand_neighbourhoods,
        ):
            if ligand_id_2 in g[ligand_id_1].neighbours():
                transform = get_transform(ligand_neighbourhood_2, ligand_neighbourhood_1)
                transforms[ligand_id_1][ligand_id_2] = transform

    return Transforms(transforms=transforms)
