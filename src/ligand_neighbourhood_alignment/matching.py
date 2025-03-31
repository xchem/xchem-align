from ligand_neighbourhood_alignment.data import Atom, AtomID, CanonicalSite, LigandNeighbourhood
from ligand_neighbourhood_alignment import dt


def match_atom(
    canonical_site_atom: Atom,
    ligand_neighbourhood_atom: Atom,
    ignore_chain=False,
) -> bool:
    id_1 = canonical_site_atom.atom_id
    id_2 = ligand_neighbourhood_atom.atom_id

    if ignore_chain:
        if id_1.atom == id_2.atom:
            if id_1.residue == id_2.residue:
                return True
    else:
        if id_1.atom == id_2.atom:
            if id_1.residue == id_2.residue:
                if id_1.chain == id_2.chain:
                    return True

    return False


def _match_atom(
    canonical_site_atom_id: tuple[str, str, str],
    ligand_neighbourhood_atom_id: tuple[str, str, str],
    ignore_chain=False,
) -> bool:
    # id_1 = canonical_site_atom.atom_id
    # id_2 = ligand_neighbourhood_atom.atom_id

    if ignore_chain:
        if canonical_site_atom_id[2] == ligand_neighbourhood_atom_id[2]:
            if canonical_site_atom_id[1] == ligand_neighbourhood_atom_id[1]:
                return True
    else:
        if canonical_site_atom_id[2] == ligand_neighbourhood_atom_id[2]:
            if canonical_site_atom_id[1] == ligand_neighbourhood_atom_id[1]:
                if canonical_site_atom_id[0] == ligand_neighbourhood_atom_id[0]:
                    return True

    return False


def match_atoms(
    site_1_atoms: dict[AtomID, Atom],
    site_2_atoms: dict[AtomID, Atom],
    min_alignable_atoms: int = 5,
) -> bool:
    # Check if there is an alignable number of atoms shared between the
    num_alignable_atoms: int = 0
    for canonical_site_atom_id, canonical_site_atom in site_1_atoms.items():
        for (
            ligand_neighbourhood_atom_id,
            ligand_neighbourhood_atom,
        ) in site_2_atoms.items():
            if match_atom(
                canonical_site_atom,
                ligand_neighbourhood_atom,
            ):
                num_alignable_atoms += 1

    if num_alignable_atoms > min_alignable_atoms:
        return True
    else:
        return False


def match_neighbourhood_to_site(
    canonical_site: CanonicalSite,
    ligand_neighbourhood: LigandNeighbourhood,
    min_alignable_atoms: int = 5,
) -> bool:
    # Check if there is an alignable number of atoms shared between the
    return match_atoms(canonical_site.atoms, ligand_neighbourhood.atoms, min_alignable_atoms)


def match_neighbourhood_to_sites(
    canonical_sites: dict[int, CanonicalSite],
    ligand_neighbourhood: LigandNeighbourhood,
) -> int | None:
    for canonical_site_id, canonical_site in canonical_sites.items():
        match: bool = match_neighbourhood_to_site(canonical_site, ligand_neighbourhood)
        if match:
            return match

    return None
