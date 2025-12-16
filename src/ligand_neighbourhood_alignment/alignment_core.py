import gemmi
from loguru import logger

from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment import structure
from ligand_neighbourhood_alignment import alignment_landmarks 
from ligand_neighbourhood_alignment.alignment_landmarks import icode_to_string


def match_atom(
    canonical_site_atom: dt.Atom,
    ligand_neighbourhood_atom: dt.Atom,
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
    site_1_atoms: dict[dt.AtomID, dt.Atom],
    site_2_atoms: dict[dt.AtomID, dt.Atom],
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
    canonical_site: dt.CanonicalSite,
    ligand_neighbourhood: dt.LigandNeighbourhood,
    min_alignable_atoms: int = 5,
) -> bool:
    # Check if there is an alignable number of atoms shared between the
    return match_atoms(canonical_site.atoms, ligand_neighbourhood.atoms, min_alignable_atoms)


def match_neighbourhood_to_sites(
    canonical_sites: dict[int, dt.CanonicalSite],
    ligand_neighbourhood: dt.LigandNeighbourhood,
) -> int | None:
    for canonical_site_id, canonical_site in canonical_sites.items():
        match: bool = match_neighbourhood_to_site(canonical_site, ligand_neighbourhood)
        if match:
            return match

    return None

def _match_cas(
    ligand_1_neighbourhood: dt.Neighbourhood,
    ligand_2_neighbourhood: dt.Neighbourhood,
    min_alignable_atoms: int = 9,  # 10 splits A71, but some things almost identical end up in different clusters
    max_alignable_rmsd: float = 1.0,
):
    alignable_cas = []
    for ligand_1_atom_id, ligand_1_atom in ligand_1_neighbourhood.atoms.items():
        for ligand_2_atom_id, ligand_2_atom in ligand_2_neighbourhood.atoms.items():
            if ligand_1_atom_id[2] == "CA":
                if _match_atom(ligand_1_atom_id, ligand_2_atom_id, ignore_chain=True):
                    alignable_cas.append(
                        (
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
                    )

    if len(alignable_cas) >= min(
        [min_alignable_atoms, len(ligand_1_neighbourhood.atoms), len(ligand_2_neighbourhood.atoms)]
    ):
        sup = gemmi.superpose_positions(
            [alignable_ca[0] for alignable_ca in alignable_cas],
            [alignable_ca[1] for alignable_ca in alignable_cas],
        )
        transform = sup.transform
        inverse_transform = sup.transform.inverse()

        rmsd = sup.rmsd
        if rmsd < max_alignable_rmsd:
            return (
                True,
                dt.Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist()),
                dt.Transform(vec=inverse_transform.vec.tolist(), mat=inverse_transform.mat.tolist()),
            )
        else:
            return False, None, None
    else:
        return False, None, None


def calculate_insertion_matching_from_landmarks(
        ref: dt.StructureLandmarks, 
        mov: dt.StructureLandmarks, 
        ref_chain: str, 
        mov_chain: str,
        debug=False
        ):
    """
    This could be a lot more efficient (and readable) it is used tables indexed to the alignment 
    space with columns for the different index types. It might need to be...
    """
    # Get the residue sequences for alignment
    ref_seq = {
        atom_id[1][1]: gemmi.one_letter_code([atom_id[1][0],]) 
        for atom_id in ref 
        if atom_id[0] == ref_chain
        }
    mov_seq = {
        atom_id[1][1]: gemmi.one_letter_code([atom_id[1][0],]) 
        for atom_id in mov 
        if atom_id[0] == mov_chain
        }
    
    # Perform the alignment
    blosum62 = gemmi.AlignmentScoring('b')
    # AA = gemmi.ResidueKind.AA
    ref_seq_sorted = [ref_seq[x] for x in ref_seq]
    mov_seq_sorted = [mov_seq[x] for x in mov_seq]
    result = gemmi.align_string_sequences(
        ref_seq_sorted,
        mov_seq_sorted,
        [],
        blosum62,
    )

    # Get native index to aligned index for the ref seq
    ref_seq_expanded = result.add_gaps(str(ref_seq_sorted), 1)
    count = 0
    ref_to_aligned = {}
    for j, item in enumerate(ref_seq_expanded):
        if item != '-':
            ref_to_aligned[count] = j
            count += 1
    
    # Get the aligned index to native index for the mov seq
    mov_seq_expanded = result.add_gaps(str(mov_seq_sorted), 2)
    mov_seq_index_to_key = {_j: _key for _j, _key in enumerate(mov_seq)}
    aligned_to_mov = {}
    count = 0
    for j, item in enumerate(mov_seq_expanded):
        if item != '-':
            aligned_to_mov[j] = mov_seq_index_to_key[count]
            count += 1
    
    # Generate the residue mapping
    count = 0
    residue_mapping = {}

    try:
        for j, ref_index in enumerate(ref_seq):
            # Get the alignment space index
            alignment_index = ref_to_aligned[j]

            # Get the corresponding mov res at that alignment index
            mov_seq_res = mov_seq_expanded[alignment_index]
            
            # Only consider matched residues
            if mov_seq_res != '-':
                # Get the mov seq index
                mov_index = aligned_to_mov[j]

                # Update the residue mapping
                residue_mapping[ref_index] = mov_index

    except Exception as e:
        print('# j, ref index')
        print([j, ref_index])
    
        print('# ref seq')
        print(ref_seq)

        print('# ref to aligned')
        print(ref_to_aligned)

        print('# mov seq expanded')
        print(mov_seq_expanded)


        raise e    

    return residue_mapping


def _get_transform_from_residues(rs: list[tuple[str, str]], srs, ssrs, other_rs=None):
    # Transform from ssrs to srs
    acs = []
    for j, resid in enumerate(rs):
        chain, num = resid
        try:
            srsr = srs[0][chain][f"{num}"][0]
            if other_rs:
                other_chain, other_num = other_rs[j]
            else:
                other_chain, other_num = chain, num
            ssrsr = ssrs[0][other_chain][f"{other_num}"][0]
            srsca = srsr["CA"][0]
            ssrsca = ssrsr["CA"][0]
            acs.append((srsca, ssrsca))
        except Exception as e:
            print(f"{chain} : {num} : {e}")

            continue

    logger.debug(f"{len(acs)}")
    if len(acs) < 3:
        print("####### SRS")

        for model in srs:
            for c in model:
                for r in c:
                    print(f"{c.name}: {r.seqid.num}")

        print("####### SSRS")
        for model in ssrs:
            for c in model:
                for r in c:
                    print(f"{c.name}: {r.seqid.num}")

        print("####### RESIDS")
        for rid in rs:
            print(f"{rid[0]} {rid[1]}")

        raise Exception()

    sup = gemmi.superpose_positions([x[0].pos for x in acs], [x[1].pos for x in acs])

    return sup.transform

def align_on_residues(
        residue_mapping,   # From reference to moving
        reference_structure, 
        moving_structure,
        ):
    """
    Take a residue mapping and generate a transform object
    """
    ref_poss = []
    mov_poss = []
    for ref_resid, mov_resid in residue_mapping.items():
        ref_res = structure.get_res_from_resid(ref_resid, reference_structure)
        ref_ca_pos = ref_res['CA'][0].pos
        mov_res = structure.get_res_from_resid(mov_resid, moving_structure)
        mov_ca_pos = mov_res['CA'][0].pos
        ref_poss.append(ref_ca_pos)
        mov_poss.append(mov_ca_pos)

    transform = gemmi.superpose_positions(ref_poss, mov_poss).transform

    return dt.Transform(
        transform.vec.tolist(),
        transform.mat.tolist(),
    )


def get_residue_mapping(  # ref-to-mov mapping
        reference_chain,
        moving_chain,
        reference_structure,
        moving_structure,
    ) -> dict[tuple[str, str]]:
    # Get landmarks from structures
    ref_landmarks = alignment_landmarks.structure_to_landmarks(reference_structure)
    mov_landmarks = alignment_landmarks.structure_to_landmarks(moving_structure)

    # Get the insertion ID mapping
    insertion_mapping = calculate_insertion_matching_from_landmarks(
        ref_landmarks,
        mov_landmarks,
        reference_chain,
        moving_chain,
    )

    # Get the resid mapping
    residue_mapping = {
        (reference_chain, str(ref_seqid)): (moving_chain, str(mov_seqid))
        for ref_seqid, mov_seqid
        in insertion_mapping.items()
        # for residue in chain
        # for chain in model
        # for model in reference_structure
        # if (chain.name == reference_chain) & (residue.seqid.num in insertion_mapping)
    }

    return residue_mapping
