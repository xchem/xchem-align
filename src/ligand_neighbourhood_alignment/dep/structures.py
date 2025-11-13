import gemmi
from loguru import logger

from ligand_neighbourhood_alignment.data import LigandNeighbourhood, ResidueID, SystemData, XtalForm
from src.ligand_neighbourhood_alignment.alignment_core import match_atom



def _get_transform():
    ...


def get_transform_from_residues(rs: list[ResidueID], srs, ssrs):
    # Transform from ssrs to srs
    acs = []
    for resid in rs:
        chain, num = resid.chain, resid.residue
        try:
            srsr = srs[0][chain][f"{num}"][0]
            ssrsr = ssrs[0][chain][f"{num}"][0]

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
            print(f"{rid.chain} {rid.residue}")

        raise Exception()

    sup = gemmi.superpose_positions([x[0].pos for x in acs], [x[1].pos for x in acs])

    return sup.transform





def get_transforms(
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

    return sup.transform, list(alignable_cas.keys())


def generate_assembly(xtalform: XtalForm, structure):
    full_st = structure.clone()
    chains_to_delete = []
    for model in full_st:
        for chain in model:
            chains_to_delete.append((str(model.num), chain.name))

    for model_name, chain_name in chains_to_delete:
        del full_st[int(model_name)-1][chain_name]

    for assembly_id, assembly in xtalform.assemblies.items():
        for generator_id, generator in assembly.generators.items():
            op = gemmi.Op(generator.triplet)
            chain_clone = structure[0][generator.chain].clone()
            for residue in chain_clone:
                for atom in residue:
                    atom_frac = structure.cell.fractionalize(atom.pos)
                    new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                    new_pos_orth = structure.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))

                    atom.pos = gemmi.Position(*new_pos_orth)
            chain_clone.name = generator.reference_chain
            full_st[0].add_chain(chain_clone)

    num_chains = 0
    for model in full_st:
        for chain in model:
            num_chains += 1
    logger.debug(f"Generated {num_chains} assembly chains")

    return full_st


def remove_non_contact_chains(assembly, neighbourhood: LigandNeighbourhood):
    chains_list = []
    for model in assembly:
        for chain in model:
            chains_list.append((str(model.num), chain.name))

    chains = list(set(chains_list))

    contact_chains_list = []
    for model in assembly:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    for atom_id, art_atom in zip(
                        neighbourhood.artefact_atom_ids,
                        neighbourhood.artefact_atoms,
                    ):
                        if atom.pos.dist(gemmi.Position(art_atom.x, art_atom.y, art_atom.z)) < 0.1:
                            contact_chains_list.append((str(model.num), chain.name))

    contact_chains = list(set(contact_chains_list))
    logger.debug(f"Num contact chains: {len(contact_chains)}")

    for model_name, chain_name in chains:
        if (model_name, chain_name) not in contact_chains:
            del assembly[int(model_name)-1][chain_name]
