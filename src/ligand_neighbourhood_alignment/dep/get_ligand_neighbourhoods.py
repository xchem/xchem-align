import gemmi
from loguru import logger

from ligand_neighbourhood_alignment import dt


# import loguru
from ligand_neighbourhood_alignment.data import (
    AssignedXtalForms,
    Atom,
    AtomID,
    Dataset,
    LigandID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Structure,
    SystemData,
    Transform,
    XtalForm,
    XtalForms,
)
from ligand_neighbourhood_alignment.structures import generate_assembly


def get_structure_fragments(dataset: Dataset, structure: Structure) -> dict[LigandID, gemmi.Residue]:
    fragments: dict[LigandID, gemmi.Residue] = {}
    # lig_number: int = 0
    for model in structure:
        for chain in model:
            for residue in chain.get_ligands():
                for lbe in dataset.ligand_binding_events.ligand_binding_events:
                    if (lbe.chain == chain.name) & (lbe.residue == residue.seqid.num):
                        ligand_id: LigandID = LigandID(
                            dtag=dataset.dtag,
                            chain=chain.name,
                            residue=lbe.residue,
                        )
                        fragments[ligand_id] = residue

    return fragments





def __get_model_and_artefact_atoms(
    residue_neighbours: dict[tuple[float, float, float], gemmi.NeighborSearch.Mark],
    structure: Structure,
) -> tuple[dict[gemmi.NeighborSearch.Mark, gemmi.CRA], dict[gemmi.NeighborSearch.Mark, gemmi.CRA]]:
    # Check each mark for its image and partition them on this
    model_atoms: dict[gemmi.NeighborSearch.Mark, gemmi.CRA] = {}
    artefact_atoms: dict[gemmi.NeighborSearch.Mark, gemmi.CRA] = {}
    for pos, cra in residue_neighbours.items():
        # Image 0 is the identity i.e. part of the normal model
        pos_gemmi = gemmi.Position(*pos)

        logger.debug(f"{cra}")
        logger.debug(f"{pos_gemmi.dist(cra.atom.pos)}")
        logger.debug(f"{pos_gemmi}")
        logger.debug(f"{cra.atom.pos}")
        if pos_gemmi.dist(cra.atom.pos) > 0.1:
            artefact_atoms[pos] = cra
        else:
            model_atoms[pos] = cra

    return model_atoms, artefact_atoms





def get_ligand_neighbourhood(
    structure: Structure,
    ns: gemmi.NeighborSearch,
    fragment: gemmi.Residue,
    min_dist: float = 0.01,
    max_dist: float = 5.0,
) -> LigandNeighbourhood:
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
            atom_id: AtomID = AtomID(
                chain=cra.chain.name,
                residue=cra.residue.seqid.num,
                atom=cra.atom.name,
            )

            # Nearest image of canon atom
            nearest_image = structure.cell.find_nearest_pbc_image(atom.pos, cra.atom.pos, neighbour.image_idx)

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
    model_atoms: dict[AtomID, Atom] = {}
    for pos, cra in _model_atoms:
        model_atom_id: AtomID = AtomID(
            chain=cra.chain.name,
            residue=cra.residue.seqid.num,
            atom=cra.atom.name,
        )
        image_transform = atom_images[model_atom_id]
        transform = Transform(
            vec=image_transform.vec.tolist(),
            mat=image_transform.mat.tolist(),
        )
        model_atoms[model_atom_id] = Atom(
            element=atom.element.name,
            atom_id=model_atom_id,
            x=pos.x,
            y=pos.y,
            z=pos.z,
            image=transform,
        )

    # Artefact atoms
    artefact_atoms: dict[AtomID, Atom] = {}
    for pos, cra in _artefact_atoms:
        artefact_atom_id: AtomID = AtomID(
            chain=cra.chain.name,
            residue=cra.residue.seqid.num,
            atom=cra.atom.name,
        )
        image_transform = atom_images[artefact_atom_id]
        transform = Transform(
            vec=image_transform.vec.tolist(),
            mat=image_transform.mat.tolist(),
        )
        artefact_atoms[artefact_atom_id] = Atom(
            element=atom.element.name,
            atom_id=artefact_atom_id,
            x=pos.x,
            y=pos.y,
            z=pos.z,
            image=transform,
        )

    # Cosntruct the neighbourhood
    ligand_neighbourhood: LigandNeighbourhood = LigandNeighbourhood(
        atom_ids=list(model_atoms.keys()),
        atoms=list(model_atoms.values()),
        artefact_atom_ids=list(artefact_atoms.keys()),
        artefact_atoms=list(artefact_atoms.values()),
    )

    return ligand_neighbourhood




def get_dataset_neighbourhoods(
    dataset: Dataset, xtalform: XtalForm, max_radius: float = 7.0
) -> dict[LigandID, LigandNeighbourhood]:
    # Load the structure
    logger.debug(dataset.pdb)
    structure: Structure = gemmi.read_structure(dataset.pdb)
    logger.debug(f"{structure.cell}")

    # Get the rest of the assembly
    assembly = generate_assembly(xtalform, structure)

    # Get the bound fragments
    fragments: dict[LigandID, gemmi.Residue] = get_structure_fragments(dataset, assembly)
    logger.debug(f"Get {len(fragments)} fragment neighbourhoods")
    logger.debug(fragments)

    # Construct the neighbourhood search
    ns: gemmi.NeighborSearch = gemmi.NeighborSearch(assembly[0], assembly.cell, max_radius).populate()

    # For each bound fragment, identify the neighbourhood atoms and
    # partition them into model and artefact
    fragment_neighbourhoods: dict[LigandID, LigandNeighbourhood] = {}
    for ligand_id, fragment in fragments.items():
        fragment_neighbourhoods[ligand_id] = get_ligand_neighbourhood(assembly, ns, fragment, max_dist=max_radius)

    return fragment_neighbourhoods


def get_ligand_neighbourhoods(
    system_data: SystemData,
    xtalforms: XtalForms,
    assigned_xtalforms: AssignedXtalForms,
) -> LigandNeighbourhoods:
    # Iterate over data, loading in structures, getting ligands for each
    # structure and finding their neighbourhoods
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood] = {}
    for dataset_id, dataset in zip(system_data.dataset_ids, system_data.datasets):
        assigned_xtalform = assigned_xtalforms.get_xtalform_id(dataset_id)
        xtalform: XtalForm = xtalforms.get_xtalform(assigned_xtalform)
        dataset_ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood] = get_dataset_neighbourhoods(
            dataset, xtalform
        )
        ligand_neighbourhoods.update(dataset_ligand_neighbourhoods)

    return LigandNeighbourhoods(
        ligand_ids=list(ligand_neighbourhoods.keys()),
        ligand_neighbourhoods=list(ligand_neighbourhoods.values()),
    )
