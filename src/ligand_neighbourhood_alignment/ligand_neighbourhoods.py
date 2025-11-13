import gemmi
from loguru import logger
from rich import print as rprint 
import numpy as np

from ligand_neighbourhood_alignment import dt
from src.ligand_neighbourhood_alignment.alignment_core import _match_cas



def _get_centroid_res(
    residues: list[tuple[str, str]],
    reference_neighbourhood: dt.Neighbourhood,
):
    res_cas = {}
    for _residue_id in residues:
        for _atom_id, _atom in reference_neighbourhood.atoms.items():
            if (_atom_id[0] == _residue_id[0]) & (_atom_id[1] == _residue_id[1]) & (_atom_id[2] == "CA"):
                res_cas[_atom_id] = _atom
    id_arr = [_atom_id for _atom_id in res_cas]
    arr = np.array([[_atom.x, _atom.y, _atom.z] for _atom in res_cas.values()])
    centroid = np.mean(arr, axis=0)
    closest = np.argmin(np.linalg.norm(arr - centroid, axis=1))
    closest_atom_id = id_arr[closest]

    return (closest_atom_id[0], closest_atom_id[1])



def _generate_assembly(xtalform: dt.XtalForm, structure, assemblies: dict[str, dt.Assembly], pdb, dataset):
    full_st = structure.clone()
    chains_to_delete = []
    for model in full_st:
        for chain in model:
            chains_to_delete.append((str(model.num), chain.name))

    for model_name, chain_name in chains_to_delete:
        del full_st[model_name][chain_name]

    cloned_chains = []
    for xtalform_assembly_id, xtalform_assembly in xtalform.assemblies.items():
        assembly = assemblies[xtalform_assembly.assembly]
        # chains = xtalform_assembly.chains
        # reference = assembly.reference
        for _biogen, _chain, _transform in zip(
            assembly.generators,
            xtalform_assembly.chains,
            xtalform_assembly.transforms,
        ):
            cloned_chains.append(_chain)

            # for generator in assembly.generators:
            #     op = gemmi.Op(generator.triplet)
            op = gemmi.Op(_transform)
            # chain_clone = structure[0][generator.chain].clone()
            try:
                chain_clone = structure[0][_chain].clone()
            except Exception as e:
                raise Exception(
                    f"An Exception occurred in generating the biological assemblies for\n"
                    f"{pdb}\n"
                    f"Based on the assembly, the expected chains were: {xtalform_assembly.chains}\n"
                    f"However the chains in the structure were: {[_x.name for _x in structure[0]]}\n"
                    "XCA does not currently handle datasets with a mis-match between the xtalform chains.\n"
                    "You should ensure that the chain names are consistent with the reference dataset for the xtalforms."
                )

            for residue in chain_clone:
                for atom in residue:
                    atom_frac = structure.cell.fractionalize(atom.pos)
                    new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                    new_pos_orth = structure.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                    atom.pos = gemmi.Position(*new_pos_orth)
            chain_clone.name = f"{_chain}~{_biogen.biomol}~{_transform}"
            full_st[0].add_chain(chain_clone)

    # Catch any ligand only chains
    for lbe in dataset.ligand_binding_events:
        _chain = lbe[1]
        if _chain not in cloned_chains:
            op = gemmi.Op("x,y,z")
            chain_clone = structure[0][_chain].clone()
            for residue in chain_clone:
                for atom in residue:
                    atom_frac = structure.cell.fractionalize(atom.pos)
                    new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                    new_pos_orth = structure.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                    atom.pos = gemmi.Position(*new_pos_orth)
            chain_clone.name = f"{_chain}~{_chain}~x,y,z"
            full_st[0].add_chain(chain_clone)
            cloned_chains.append(_chain)

    chains = []
    num_chains = 0
    for model in full_st:
        for chain in model:
            num_chains += 1
            chains.append(chain.name)
    logger.debug(f"Generated {num_chains} assembly chains")
    logger.debug(f"Chain names are: {[x.name for x in full_st[0]]}")
    # print(chains)

    return full_st


def _get_structure_fragments(dataset: dt.Dataset, structure, version):
    fragments: dict[dt.LigandNeighbourhoodID, gemmi.Residue] = {}
    for model in structure:
        for chain in model:
            source_chain, biomol_chain, transform = chain.name.split("~")
            for residue in chain:  # .get_ligands():
                for lbe in dataset.ligand_binding_events:
                    if (
                        (str(lbe[2]) == str(residue.seqid.num))
                        & (str(lbe[1]) == str(source_chain))
                        & (transform == "x,y,z")
                    ):
                        ligand_id = (dataset.dtag, str(lbe[1]), str(lbe[2]), str(lbe[3]), str(version))
                        fragments[ligand_id] = residue

    return fragments




def get_model_and_artefact_atoms(
    residue_neighbours: list[tuple[gemmi.Position, gemmi.CRA]],
    structure: dt.Structure,
    fragment,
) -> tuple[list[tuple[gemmi.Position, gemmi.CRA]], list[tuple[gemmi.Position, gemmi.CRA]]]:
    # Check each mark for its image and partition them on this
    model_atoms: list[tuple[gemmi.Position, gemmi.CRA]] = []
    artefact_atoms: list[tuple[gemmi.Position, gemmi.CRA]] = []
    possible_artefact_atoms = []
    for pos, cra in residue_neighbours:
        # Image 0 is the identity i.e. part of the normal model
        canon_atom = cra.atom
        # Possible artefact: Could be sym, ncs or translation
        if canon_atom.pos.dist(pos) > 0.1:
            possible_artefact_atoms.append((pos, cra))

        # Canonical/model atom confirnmed: just see if already handled
        else:
            # Check there is not already a nearby atom in there
            if all([model_atom[0].dist(pos) > 0.1 for model_atom in model_atoms]):
                model_atoms.append((pos, cra))

    for pos, cra in possible_artefact_atoms:
        # Check it isn't a ncs image by seeing if it overlays a model atom
        if all([model_atom[0].dist(pos) > 0.1 for model_atom in model_atoms]):
            if all([model_atom[0].dist(pos) > 0.1 for model_atom in artefact_atoms]):
                artefact_atoms.append((pos, cra))

    updated_model_atoms = [
        atom for atom in model_atoms if all([fragment_atom.pos.dist(atom[0]) > 0.1 for fragment_atom in fragment])
    ]

    updated_artefact_atoms = [
        atom for atom in artefact_atoms if all([fragment_atom.pos.dist(atom[0]) > 0.1 for fragment_atom in fragment])
    ]

    # return model_atoms, artefact_atoms
    return updated_model_atoms, updated_artefact_atoms

def _get_ligand_neighbourhood(
    structure,
    ns: gemmi.NeighborSearch,
    fragment: gemmi.Residue,
    min_dist: float = 0.01,
    max_dist: float = 5.0,
):
    # For each atom, get the neighbouring atoms, and filter them on their
    # real space position
    residue_neighbours: list[tuple[gemmi.Position, gemmi.CRA]] = []

    atom_images = {}
    for atom in fragment:
        # Get the atom neighbours (as arbitrary image marks)
        atom_neighbours: list[gemmi.NeighborSearch.Mark] = ns.find_neighbors(
            atom,
            min_dist=min_dist,
            max_dist=max_dist,
        )
        for neighbour in atom_neighbours:
            # Get positions of the marks nearest image to canon atom
            cra = neighbour.to_cra(structure[0])
            atom_id: tuple[str, str, str] = (
                str(cra.chain.name),
                str(cra.residue.seqid.num),
                str(cra.atom.name),
            )

            # Nearest image of canon atom
            nearest_image = structure.cell.find_nearest_pbc_image(
                atom.pos,
                cra.atom.pos,
                neighbour.image_idx,
            )

            # Get canon atom pos as tractional
            fpos = structure.cell.fractionalize(cra.atom.pos)

            # Get transform that generates image from canon
            ftransform = ns.get_image_transformation(neighbour.image_idx)
            atom_images[atom_id] = ftransform

            # Apply the canon -> image and image -> pbc image transforms
            fpos_trans = ftransform.apply(fpos)
            fpos_trans.x = fpos_trans.x + nearest_image.pbc_shift[0]
            fpos_trans.y = fpos_trans.y + nearest_image.pbc_shift[1]
            fpos_trans.z = fpos_trans.z + nearest_image.pbc_shift[2]

            pos = structure.cell.orthogonalize(fpos_trans)

            residue_neighbours.append((pos, cra))

    # # Seperate out model and artefact atoms
    _model_atoms, _artefact_atoms = get_model_and_artefact_atoms(residue_neighbours, structure, fragment)
    logger.debug(f"Got {len(_model_atoms)} model atoms")
    logger.debug(f"Got {len(_artefact_atoms)} artefact atoms")

    # Model atoms
    model_atoms: dict[tuple[str, str, str], dt.Atom] = {}
    for pos, cra in _model_atoms:
        # cra = atom.to_cra(structure[0])
        model_atom_id: tuple[str, str, str] = (
            str(cra.chain.name),
            str(cra.residue.seqid.num),
            str(cra.atom.name),
        )
        image_transform = atom_images[model_atom_id]
        model_atom_id: tuple[str, str, str] = (
            str(cra.chain.name).split("~", maxsplit=1)[0],
            str(cra.residue.seqid.num),
            str(cra.atom.name),
        )
        transform = dt.Transform(
            vec=image_transform.vec.tolist(),
            mat=image_transform.mat.tolist(),
        )
        model_atoms[model_atom_id] = dt.Atom(
            element=cra.atom.element.name,
            x=pos.x,
            y=pos.y,
            z=pos.z,
            image=transform,
        )

    # Artefact atoms
    artefact_atoms: dict[tuple[str, str, str], dt.Atom] = {}
    for pos, cra in _artefact_atoms:
        artefact_atom_id: tuple[str, str, str] = (
            str(cra.chain.name),
            str(cra.residue.seqid.num),
            str(cra.atom.name),
        )
        image_transform = atom_images[artefact_atom_id]
        artefact_atom_id: tuple[str, str, str] = (
            str(cra.chain.name).split("~", maxsplit=1)[0],
            str(cra.residue.seqid.num),
            str(cra.atom.name),
        )
        transform = dt.Transform(
            vec=image_transform.vec.tolist(),
            mat=image_transform.mat.tolist(),
        )
        artefact_atoms[artefact_atom_id] = dt.Atom(
            element=cra.atom.element.name,
            x=pos.x,
            y=pos.y,
            z=pos.z,
            image=transform,
        )

    # Cosntruct the neighbourhood
    ligand_neighbourhood = dt.Neighbourhood(model_atoms, artefact_atoms)

    return ligand_neighbourhood


def _get_dataset_neighbourhoods(
    dataset: dt.Dataset, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly], version, max_radius: float = 9.0
) -> dict[dt.LigandNeighbourhoodID, dt.Neighbourhood]:
    # Load the structure
    logger.debug(dataset.pdb)
    structure = gemmi.read_structure(dataset.pdb)
    logger.debug(f"{structure.cell}")

    # Get the rest of the assembly
    assembly = _generate_assembly(xtalform, structure, assemblies, dataset.pdb, dataset)

    # Get the bound fragments
    fragments: dict[dt.LigandNeighbourhoodID, gemmi.Residue] = _get_structure_fragments(dataset, assembly, version)
    logger.debug(f"Get {len(fragments)} fragment neighbourhoods")
    logger.debug(fragments)

    # Construct the neighbourhood search
    ns: gemmi.NeighborSearch = gemmi.NeighborSearch(
        assembly[0],
        assembly.cell,
        max_radius,
    ).populate()

    # For each bound fragment, identify the neighbourhood atoms and
    # partition them into model and artefact
    fragment_neighbourhoods: dict[dt.LigandNeighbourhoodID, dt.Neighbourhood] = {}
    for ligand_id, fragment in fragments.items():
        fragment_neighbourhoods[ligand_id] = _get_ligand_neighbourhood(
            assembly,
            ns,
            fragment,
            max_dist=max_radius,
        )

    return fragment_neighbourhoods


def _get_neighbourhoods(
    dataset: dt.Dataset,
    xtalform: dt.XtalForm,
    assemblies: dict[str, dt.Assembly],
    version,
):
    dataset_ligand_neighbourhoods: dict[dt.LigandNeighbourhoodID, dt.Neighbourhood] = _get_dataset_neighbourhoods(
        dataset, xtalform, assemblies, version
    )
    return dataset_ligand_neighbourhoods

def _update_ligand_neighbourhood_transforms(
    ligand_neighbourhood_transforms: dict[tuple[dt.LigandNeighbourhoodID, dt.LigandNeighbourhoodID], dt.Transform],
    lid: dt.LigandNeighbourhoodID,
    ligand_neighbourhoods: dict[dt.LigandNeighbourhoodID, dt.Neighbourhood],
    structures,
):
    # For each other ligand neighbourhood, see if atoms match and then if so record the transform that relates
    # them and its inverse

    ligand_1_id = lid
    ligand_1_neighbourhood = ligand_neighbourhoods[lid]
    matches = []
    for (
        ligand_2_id,
        ligand_2_neighbourhood,
    ) in ligand_neighbourhoods.items():
        # See if atoms match - transform is frame 2 to frame 1
        ca_match, transform, inverse_transform = _match_cas(ligand_1_neighbourhood, ligand_2_neighbourhood)

        if ca_match:
            ligand_neighbourhood_transforms[(ligand_1_id, ligand_2_id)] = transform
            ligand_neighbourhood_transforms[(ligand_2_id, ligand_1_id)] = inverse_transform
            matches.append(ligand_2_id)

    if len(matches) == 0:
        rprint(f"No Matches For {ligand_1_id}! No alignments will be generated!")
