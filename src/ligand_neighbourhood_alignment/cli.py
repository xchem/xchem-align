import json
import os
import sys

# import os
import subprocess

# import sys
from pathlib import Path

import fire
import numpy as np
import pandas as pd
import yaml
from loguru import logger

import networkx as nx

from rich import print

# in the original code, the imports were like this:
# from rich import print as rprint

# logger.remove()  # for someone not familiar with the lib, whats going on here?
# logger.add(sys.stdout, level="WARNING")
# from rich import print
# I wonder if this has some significance?

import gemmi

from ligand_neighbourhood_alignment import dt

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.align_xmaps import _align_xmaps

from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import _get_ligand_neighbourhood


# from ligand_neighbourhood_alignment.get_system_sites import get_system_sites
from ligand_neighbourhood_alignment.build_alignment_graph import build_alignment_graph
from ligand_neighbourhood_alignment.data import (  # save_xtalforms,
    Assemblies,
    AssignedXtalForms,
    CanonicalSites,
    ChainOutput,
    ConformerSites,
    Dataset,
    DatasetID,
    DatasetOutput,
    Datasource,
    LigandID,
    LigandNeighbourhoods,
    LigandOutput,
    Options,
    Output,
    PanDDA,
    SiteTransforms,
    SystemData,
    Transforms,
    XtalForms,
    read_assigned_xtalforms,
    read_canonical_sites,
    read_graph,
    read_neighbourhoods,
    read_output,
    read_site_transforms,
    read_structures,
    read_system_data,
    read_transforms,
    read_xtalforms,
    save_assigned_xtalforms,
    save_canonical_sites,
    save_data,
    save_output,
)
from ligand_neighbourhood_alignment.generate_aligned_structures import _align_structures_from_sites
from ligand_neighbourhood_alignment.generate_sites_from_components import (  # get_xtalform_sites_from_canonical_sites,
    _generate_sites_from_components,
    get_components,
    get_conformer_sites_from_components,
    get_site_transforms,
    get_sites_from_conformer_sites,
    get_structures,
    get_subsite_transforms,
    _update_conformer_site_transforms,
    _update_canonical_site_transforms,
    _update_reference_structure_transforms,
)
from ligand_neighbourhood_alignment.get_alignability import get_alignability, _update_ligand_neighbourhood_transforms
from ligand_neighbourhood_alignment.get_graph import get_graph
from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from ligand_neighbourhood_alignment.make_data_json import (
    get_ligand_binding_events_from_panddas,
    get_ligand_binding_events_from_structure,
    make_data_json_from_pandda_dir,
)
from ligand_neighbourhood_alignment.generate_aligned_structures import _align_structure, _align_reference_structure
from ligand_neighbourhood_alignment.align_xmaps import read_xmap, read_xmap_from_mtz, __align_xmap
from ligand_neighbourhood_alignment import alignment_heirarchy


logger.remove()  # for someone not familiar with the lib, whats going on here?
logger.add(sys.stdout, level="WARNING")


def cas_ligands():
    return "\tgraphics_to_ca_plus_ligands_sec_struct_representation(p) \n"


def _change_sites_reference(_source_dir: Path, site_id: int):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    sites.reference_site_id = site_id
    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _change_site_reference(_source_dir: Path, site_id: int, subsite_id: int):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    site = sites.get_site(site_id)
    site.reference_subsite_id = subsite_id
    new_reference_subsite = site.get_subsite(subsite_id)
    site.reference_subsite = new_reference_subsite
    site.reference_ligand_id = new_reference_subsite.reference_ligand_id
    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _change_subsite_reference(
    _source_dir: Path,
    site_id: int,
    subsite_id: int,
    dtag: int,
    chain: str,
    residue: int,
):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    site = sites.get_site(site_id)
    site.reference_subsite_id = subsite_id
    subsite = site.get_subsite(subsite_id)
    new_lid = LigandID(dtag=dtag, chain=chain, residue=residue)
    if new_lid not in subsite.members:
        raise Exception(f"LigandID {new_lid} not in {subsite.members}")
    subsite.reference_ligand_id = new_lid

    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _add_model_building_dir_to_system_data(system_data: SystemData, _data_source_dir: Path):
    datasource = Datasource(path=str(_data_source_dir), datasource_type="model_building")

    if not system_data.datasources:
        logger.info("No Datasources: Creating new list!")
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
            _datasource for _datasource in system_data.datasources if _datasource.path != str(_data_source_dir)
        ] + [
            datasource,
        ]
        system_data.datasources = new_datasources

    return system_data


def _add_model_building_dir(_source_dir: Path, _data_source_dir: Path):
    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")
    system_data = read_system_data(_source_dir)

    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")

    system_data = _add_model_building_dir_to_system_data(system_data, _data_source_dir)

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are now: {datasource_paths}")


def _add_manual_dir_to_system_data(system_data: SystemData, _data_source_dir: Path):
    datasource = Datasource(path=str(_data_source_dir), datasource_type="manual")

    if not system_data.datasources:
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
            _datasource for _datasource in system_data.datasources if _datasource.path != str(_data_source_dir)
        ] + [
            datasource,
        ]
        system_data.datasources = new_datasources

    return system_data


def _add_manual_dir(_source_dir: Path, _data_source_dir: Path):
    system_data = read_system_data(_source_dir)

    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")

    system_data = _add_manual_dir_to_system_data(system_data, _data_source_dir)

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")


def _add_pandda_to_system_data(system_data: SystemData, _pandda_dir: Path):
    analyses_dir: Path = _pandda_dir / constants.PANDDA_ANALYSES_DIR
    event_table_path: Path = analyses_dir / constants.PANDDA_EVENTS_INSPECT_TABLE_PATH

    if event_table_path.exists():
        pandda = PanDDA(path=str(_pandda_dir), event_table_path=str(event_table_path))

        if not system_data.panddas:
            system_data.panddas = [
                pandda,
            ]
        else:
            new_panddas = [_pandda for _pandda in system_data.panddas if _pandda.path != str(_pandda_dir)] + [
                pandda,
            ]
            system_data.panddas = new_panddas

    else:
        raise Exception(f"No event table at: {event_table_path}")

    return system_data


def _add_pandda(_source_dir: Path, _pandda_dir: Path):
    system_data = read_system_data(_source_dir)

    system_data = _add_pandda_to_system_data(system_data, _pandda_dir)

    save_data(system_data, _source_dir)

    logger.info(f"Added PanDDA {_pandda_dir} to panddas")
    pandda_paths = [_pandda.path for _pandda in system_data.panddas]
    logger.debug(f"PanDDAs are: {pandda_paths}")


def _add_data_to_system_data(system_data):
    # Get the PanDDA event tables
    pandda_event_tables = {pandda.path: pd.read_csv(pandda.event_table_path) for pandda in system_data.panddas}
    logger.info(f"Read {len(pandda_event_tables)} PanDDA event tables")

    # Get the
    datasets = {}
    for datasource in system_data.datasources:
        logger.info(f"Parsing datasource: {datasource.path}")
        if datasource.datasource_type == "model_building":
            for model_dir in Path(datasource.path).glob("*"):
                dtag = model_dir.name
                dataset_id = DatasetID(dtag=dtag)
                if dataset_id in datasets:
                    st = f"Dataset ID {dataset_id} already found! Using new!"
                    logger.warning(st)
                    continue

                pdb = model_dir / constants.MODEL_DIR_PDB
                xmap = model_dir / constants.MODEL_DIR_XMAP
                mtz = model_dir / constants.MODEL_DIR_MTZ
                if not pdb.exists():
                    continue

                ligand_binding_events = get_ligand_binding_events_from_panddas(pandda_event_tables, pdb, dtag)
                if len(ligand_binding_events.ligand_ids) == 0:
                    logger.warning(f"Dataset {dtag} has no ligand binding events!")
                    continue
                dataset = Dataset(
                    dtag=dtag,
                    pdb=str(pdb),
                    xmap=str(xmap),
                    mtz=str(mtz),
                    ligand_binding_events=ligand_binding_events,
                )
                # dataset_ids.append(dataset_id)
                datasets[dataset_id] = dataset
                logger.debug(f"Added dataset: {dataset_id}")

        elif datasource.datasource_type == "manual":
            for model_dir in Path(datasource.path).glob("*"):
                dtag = model_dir.name
                dataset_id = DatasetID(dtag=dtag)

                if dataset_id in datasets:
                    st = f"Dataset ID {dataset_id} already found! Using new!"
                    logger.warning(st)
                try:
                    pdb = next(model_dir.glob("*.pdb"))
                except Exception:
                    raise Exception(f"Could not find pdb in dir: {model_dir}")
                try:
                    xmap = next(model_dir.glob("*.ccp4"))
                except Exception as e:
                    print(e)
                    xmap = None
                    logger.warning("No xmap!")
                try:
                    mtz = next(model_dir.glob("*.mtz"))
                except Exception as e:
                    print(e)
                    mtz = None
                    logger.warning("No mtz!")

                ligand_binding_events = get_ligand_binding_events_from_structure(pdb, xmap, dtag)
                if len(ligand_binding_events.ligand_ids) == 0:
                    logger.warning(f"Dataset {dtag} has no ligand binding events!")
                    continue
                dataset = Dataset(
                    dtag=dtag,
                    pdb=str(pdb),
                    xmap=str(xmap),
                    mtz=str(mtz),
                    ligand_binding_events=ligand_binding_events,
                )
                datasets[dataset_id] = dataset
                logger.debug(f"Added dataset: {dataset_id}")
        else:
            raise Exception(f"Source type {datasource.datasource_type} unknown!")

    system_data.dataset_ids = list(datasets.keys())
    system_data.datasets = list(datasets.values())

    return system_data


def _parse_data_sources(_source_dir: Path):
    system_data = read_system_data(_source_dir)

    system_data = _add_data_to_system_data(system_data)

    save_data(system_data, _source_dir)

    logger.info(f"Found {len(system_data.dataset_ids)} datasets!")


def save_schema(model, path):
    with open(path / model.__name__, "w") as f:
        f.write(model.schema_json(indent=2))


def get_closest_xtalform(xtalforms: XtalForms, structures, dataset_id):
    structure = structures[dataset_id]
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    for xtalform_id, xtalform in xtalforms.iter():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        if ref_spacegroup != structure_spacegroup:
            continue

        deltas = np.array(
            [
                structure_cell.a / ref_structure_cell.a,
                structure_cell.b / ref_structure_cell.b,
                structure_cell.c / ref_structure_cell.c,
                structure_cell.alpha / ref_structure_cell.alpha,
                structure_cell.beta / ref_structure_cell.beta,
                structure_cell.gamma / ref_structure_cell.gamma,
            ]
        )
        xtalform_deltas[xtalform_id] = deltas

    if len(xtalform_deltas) == 0:
        return None, None

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform]


def _get_assigned_xtalforms(system_data, xtalforms):
    structures = read_structures(system_data)

    dataset_ids = []
    xtalform_ids = []
    for dataset_id, dataset in system_data.iter():
        closest_xtalform_id, deltas = get_closest_xtalform(xtalforms, structures, dataset_id)

        if (closest_xtalform_id is None) & (deltas is None):
            logger.info(f"No reference in same spacegroup for: {dataset_id}")
            logger.info(f"Structure path is: {dataset.pdb}")
            raise Exception()

        if np.any(deltas > 1.1) | np.any(deltas < 0.9):
            logger.info(f"No reference for dataset: {dataset_id}")
            logger.info(f"Deltas to closest unit cell are: {deltas}")
            logger.info(f"Structure path is: {dataset.pdb}")

            raise Exception()

        dataset_ids.append(dataset_id)
        xtalform_ids.append(closest_xtalform_id)

    assigned_xtalforms = AssignedXtalForms(dataset_ids=dataset_ids, xtalform_ids=xtalform_ids)

    return assigned_xtalforms


def _assign_xtalforms(
    _source_dir: Path,
    assemblies: Assemblies,
    xtalforms: XtalForms,
    system_data: SystemData,
):
    assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

    save_assigned_xtalforms(_source_dir, assigned_xtalforms)

    return assigned_xtalforms


def _get_structures(datasets):
    structures = {}
    for dtag, dataset in datasets.items():
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        # if structure.cell.a == 0.0:
        cell = structure.cell

        if cell.a == 1.0:
            poss = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            pos = atom.pos
                            poss.append((pos.x, pos.y, pos.z))
            pos_array = np.array(poss)

            cell_lengths = np.max(pos_array, axis=0) - np.min(pos_array, axis=0)

            structure.cell = gemmi.UnitCell(cell_lengths[0], cell_lengths[1], cell_lengths[2], 90.0, 90.0, 90.0)
        # print(structure.cell.a)
        structures[dataset.dtag] = structure

    return structures


def _get_dataset_protein_chains(structure):
    protein_chains = []
    for model in structure:
        for chain in model:
            protein = False
            for residue in chain:
                if residue.name in constants.RESIDUE_NAMES:
                    protein = True
            if protein:
                protein_chains.append(chain.name)

    return protein_chains


def _get_closest_xtalform(xtalforms: dict[str, dt.XtalForm], structure, structures):
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    for xtalform_id, xtalform in xtalforms.items():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        # Check they are in the same spacegroup
        if ref_spacegroup != structure_spacegroup:
            continue

        # Check they have the same protein chains
        xtalform_protein_chains = [
            _chain for xtalform_assembly in xtalform.assemblies.values() for _chain in xtalform_assembly.chains
        ]
        dataset_protein_chains = _get_dataset_protein_chains(structure)
        # print(f'Xtalform chains: {set(xtalform_protein_chains)}')
        # print(f'Dataseet chains: {set(dataset_protein_chains)}')
        # print(set(dataset_protein_chains) == set(xtalform_protein_chains))

        if set(dataset_protein_chains) != set(xtalform_protein_chains):
            continue

        deltas = np.array(
            [
                structure_cell.a / ref_structure_cell.a,
                structure_cell.b / ref_structure_cell.b,
                structure_cell.c / ref_structure_cell.c,
                structure_cell.alpha / ref_structure_cell.alpha,
                structure_cell.beta / ref_structure_cell.beta,
                structure_cell.gamma / ref_structure_cell.gamma,
            ]
        )
        xtalform_deltas[xtalform_id] = deltas

    if len(xtalform_deltas) == 0:
        return None, None

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform]


def _assign_dataset(dataset, assemblies, xtalforms, structure, structures):
    closest_xtalform_id, deltas = _get_closest_xtalform(
        xtalforms,
        structure,
        structures,
    )

    if (closest_xtalform_id is None) & (deltas is None):
        logger.info(f"No reference in same spacegroup as: {dataset.dtag}")
        logger.info(f"Structure path is: {dataset.pdb}")
        raise Exception(f"No reference in same spacegroup as: {dataset.dtag}\nStructure path is: {dataset.pdb}")

    if np.any(deltas > 1.1) | np.any(deltas < 0.9):
        logger.info(f"No reference for dataset: {dataset.dtag}")
        logger.info(f"Deltas to closest unit cell are: {deltas}")
        logger.info(f"Structure path is: {dataset.pdb}")

        raise Exception(
            f"No reference for dataset: {dataset.dtag}\nDeltas to closest unit cell in {closest_xtalform_id} are: {deltas}\nStructure path is: {dataset.pdb}"
        )

    return closest_xtalform_id


def _save_assignments(fs_model: dt.FSModel, dataset_assignments: dict[str, str]):
    with open(fs_model.dataset_assignments, "w") as f:
        yaml.safe_dump(dataset_assignments, f)


def _generate_assembly(xtalform: dt.XtalForm, structure, assemblies: dict[str, dt.Assembly], pdb, dataset):
    full_st = structure.clone()
    chains_to_delete = []
    for model in full_st:
        for chain in model:
            chains_to_delete.append((model.name, chain.name))

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
    fragments: dict[tuple[str, str, str, str], gemmi.Residue] = {}
    # lig_number: int = 0
    for model in structure:
        for chain in model:
            source_chain, biomol_chain, transform = chain.name.split("~")
            for residue in chain:  # .get_ligands():
                for lbe in dataset.ligand_binding_events:
                    # if (
                    #     (residue.name == "LIG")
                    #     & (lbe.chain == chain.name)
                    #     & (lbe.residue == residue.seqid.num)
                    # ):
                    # if (lbe.chain == chain.name) & (lbe.residue == residue.seqid.num):
                    #     ligand_id = (dataset.dtag, str(chain.name), str(lbe.residue),)
                    #     fragments[ligand_id] = residue
                    # lig_number = lig_number + 1
                    if (
                        (str(lbe[2]) == str(residue.seqid.num))
                        & (str(lbe[1]) == str(source_chain))
                        & (transform == "x,y,z")
                    ):
                        ligand_id = (dataset.dtag, str(lbe[1]), str(lbe[2]), str(version))
                        fragments[ligand_id] = residue

    return fragments


def _get_dataset_neighbourhoods(
    dataset: dt.Dataset, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly], version, max_radius: float = 9.0
) -> dict[tuple[str, str, str, str], dt.Neighbourhood]:
    # Load the structure
    logger.debug(dataset.pdb)
    structure = gemmi.read_structure(dataset.pdb)
    logger.debug(f"{structure.cell}")

    # Get the rest of the assembly
    assembly = _generate_assembly(xtalform, structure, assemblies, dataset.pdb, dataset)

    # Get the bound fragments
    fragments: dict[tuple[str, str, str, str], gemmi.Residue] = _get_structure_fragments(dataset, assembly, version)
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
    fragment_neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood] = {}
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
    dataset_ligand_neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood] = _get_dataset_neighbourhoods(
        dataset, xtalform, assemblies, version
    )
    return dataset_ligand_neighbourhoods


def _save_neighbourhoods(
    fs_model: dt.FSModel,
    ligand_neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
):
    with open(fs_model.ligand_neighbourhoods, "w") as f:
        dic = {}
        for ligand_id, neighbourhood in ligand_neighbourhoods.items():
            dic["/".join(ligand_id)] = neighbourhood.to_dict()
        yaml.safe_dump(dic, f)


def _save_ligand_neighbourhood_transforms(fs_model, ligand_neighbourhood_transforms):
    with open(fs_model.ligand_neighbourhood_transforms, "w") as f:
        dic = {}
        for (to_ligand_id, from_ligand_id), transform in ligand_neighbourhood_transforms.items():
            key = "~".join(["/".join(to_ligand_id), "/".join(from_ligand_id)])
            dic[key] = transform.to_dict()
        yaml.safe_dump(dic, f)


def _update_graph(
    alignability_graph,
    ligand_neighbourhoods,
    ligand_neighbourhood_transforms,
):
    nodes = alignability_graph.nodes
    edges = alignability_graph.edges

    for ligand_id in ligand_neighbourhoods:
        if ligand_id not in nodes:
            alignability_graph.add_node(ligand_id)

    for to_ligand_id, from_ligand_id in ligand_neighbourhood_transforms:
        # if to_ligand_id not in nodes:
        #     alignability_graph.add_node(to_ligand_id)
        # if from_ligand_id not in nodes:
        #     alignability_graph.add_node(from_ligand_id)
        if (to_ligand_id, from_ligand_id) not in edges:
            alignability_graph.add_edge(to_ligand_id, from_ligand_id)


def _save_graph(fs_model, alignability_graph):
    graph_for_output = nx.relabel_nodes(alignability_graph, {x: "/".join(x) for x in alignability_graph})
    nx.write_gml(
        graph_for_output,
        str(fs_model.alignability_graph),
        # stringizer=lambda x: "/".join(x),
    )


# def _get_connected_components(alignability_graph):
#     cliques = list(nx.connected_components(alignability_graph))
#     return cliques


def _get_connected_components(alignability_graph, clusters, max_path_length=2):
    """
    Construct neighbourhoods around the most connected neighbourhoods by some max path length,


    """

    # Get the graph of short paths
    path = dict(nx.all_pairs_shortest_path(alignability_graph))
    path_lengths = {(source, target): len(path[source][target]) for source in path for target in path[source]}
    H = nx.Graph()
    for node in path:
        H.add_node(node)
    for source, target in path_lengths:
        if path_lengths[(source, target)] <= max_path_length:
            H.add_edge(source, target)

    #
    degrees = dict(nx.degree(H))

    # Replay cluster cores
    used = [member for cluster in clusters for member in cluster]
    for x in clusters:
        for target in H.nodes:
            if target in used:
                continue
            if (x, target) not in path_lengths:
                continue
            if path_lengths[(x, target)] <= 2:
                used.append(target)
                clusters[x].append(target)

    # Now go through any new ligands that are not yet connected, constructing clusters for them
    for x in sorted(degrees, key=lambda _x: degrees[_x], reverse=True):
        if x in used:
            continue
        clusters[x] = []
        # print(f"f{x} : {degrees[x]}")

        # for n in G.neighbors(x):
        # used.append(n)
        for target in H.nodes:
            if target in used:
                continue
            if (x, target) not in path_lengths:
                continue
            if path_lengths[(x, target)] <= 2:
                used.append(target)
                clusters[x].append(target)

    return clusters


def _update_conformer_sites(
    conformer_sites: dict[str, dt.ConformerSite],
    connected_component_id: tuple[str, str, str, str],
    connected_component: list[tuple[str, str, str, str]],
    neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
    structures,
    assemblies,
    xtalforms,
    xtalform_assignments,
):
    matched = False
    # Check each old conformer site for overlap in membership, and if so update its members
    for conformer_site_id, conformer_site in conformer_sites.items():
        num_overlaps = set(connected_component).intersection(set(conformer_site.members))
        if len(num_overlaps) > 0:
            matched = True
            # Match, add each new ligand id to the conformer site's members
            for lid in connected_component:
                if lid not in conformer_site.members:
                    conformer_site.members.append(lid)

    # Otherwise create a new conformer site
    if not matched:
        residues = []
        residues_aligned = []
        for lid in connected_component:
            st = structures[lid[0]]
            if lid != connected_component_id:
                continue
            for atom_id in neighbourhoods[lid].atoms:
                biochain = alignment_heirarchy._chain_to_biochain(
                    atom_id[0], xtalforms[xtalform_assignments[lid[0]]], assemblies
                )
                residues.append((atom_id[0], atom_id[1], st[0][atom_id[0]][atom_id[1]][0].name))
                residues_aligned.append((biochain, atom_id[1], st[0][atom_id[0]][atom_id[1]][0].name))
        conformer_site = dt.ConformerSite(
            [x for x in set(residues)],
            [x for x in set(residues_aligned)],
            connected_component,
            # [x for x in connected_component][0]
            connected_component_id,
        )
        conformer_site_id = "+".join(conformer_site.reference_ligand_id)
        conformer_sites[conformer_site_id] = conformer_site


def _save_connected_components(fs_model, connected_components):
    with open(fs_model.connected_components, "w") as f:
        dic = {}
        for connected_component_reference, connected_component in connected_components.items():
            dic["+".join(connected_component_reference)] = ["+".join(member) for member in connected_component]
        yaml.safe_dump(dic, f, sort_keys=False)


def _save_conformer_sites(fs_model: dt.FSModel, conformer_sites: dict[str, dt.ConformerSite]):
    with open(fs_model.conformer_sites, "w") as f:
        dic = {}
        for conformer_site_id, conformer_site in conformer_sites.items():
            dic[conformer_site_id] = conformer_site.to_dict()
        yaml.safe_dump(dic, f)


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


def _update_canonical_sites(
    canonical_sites: dict[str, dt.CanonicalSite],
    conformer_site: dt.ConformerSite,
    conformer_site_id,
    neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
    min_shared_residues=6,
):
    if len(canonical_sites) != 0:
        global_reference_dtag = [x for x in canonical_sites.values()][0].global_reference_dtag
    else:
        global_reference_dtag = conformer_site.reference_ligand_id[0]

    # If conformer site already in a canonical site skip
    for canonical_site_id, canonical_site in canonical_sites.items():
        if conformer_site_id in canonical_site.conformer_site_ids:
            return

            # Check each canonical site to see if conformer site already in it and if not
    # whether it shares enough residues to now be added
    matched = False
    conformer_site_residues = [(residue[1], residue[2]) for residue in conformer_site.residues]
    for canonical_site_id, canonical_site in canonical_sites.items():
        if matched:
            continue
        canonical_site_residues = [(residue[1], residue[2]) for residue in canonical_site.residues]
        if conformer_site_id not in canonical_site.conformer_site_ids:
            v = set(canonical_site_residues).intersection(set(conformer_site_residues))
            # if len(v) >= min(min_shared_residues, int((3/4)*len(canonical_site_residues))):
            if len(v) >= 0.75 * len(canonical_site_residues):
                # Matched!
                matched = True
                canonical_site.conformer_site_ids.append(conformer_site_id)

    # If not matched to any existing canonical site create a new one
    if not matched:
        centroid_res = _get_centroid_res(conformer_site.residues, neighbourhoods[conformer_site.reference_ligand_id])
        canonical_site = dt.CanonicalSite(
            [
                conformer_site_id,
            ],
            conformer_site.residues,
            conformer_site_id,
            global_reference_dtag,
            (
                conformer_site.reference_ligand_id[0],
                centroid_res[0],
                centroid_res[1],
                conformer_site.reference_ligand_id[1],
            ),
        )

        canonical_site_id = conformer_site_id
        canonical_sites[canonical_site_id] = canonical_site


def _save_canonical_sites(fs_model, canonical_sites: dict[str, dt.CanonicalSite]):
    with open(fs_model.canonical_sites, "w") as f:
        dic = {}
        for canonical_site_id, canonical_site in canonical_sites.items():
            dic[canonical_site_id] = canonical_site.to_dict()
        yaml.safe_dump(dic, f)


def _get_dist(pos_1, pos_2):
    return np.linalg.norm(np.array([pos_1.x, pos_1.y, pos_1.z]) - np.array([pos_2.x, pos_2.y, pos_2.z]))


def _crystalform_incremental_cluster(observation_centroid_residues, xtalform_sites, neighbourhoods, cutoff=10.0):
    # Check if each observation centroid is within 10A of an xtalform site centroid, and if so assign that dataset to it
    # If not, select the remainer with the most members within 10A and form a new "site" from it
    # Continue until all observations are assigned

    # Get CA positions
    centroid_ca_positions = {
        observation_id: neighbourhoods[observation_id].atoms[(centroid_res[0], centroid_res[1], 'CA')]
        for observation_id, centroid_res in observation_centroid_residues.items()
    }

    # Identify current xtalform centre residues (for the considered canonical site)
    centre_residues_positions = {
        xtalform_site_id: centroid_ca_positions[tuple(xtalform_site_id.split('/'))]
        for xtalform_site_id in xtalform_sites
    }

    # While observations remain to be assigned
    assignments = {
        xtalform_site_id: [_x for _x in xtalform_sites[xtalform_site_id].members]
        for xtalform_site_id in xtalform_sites
    }
    assigned_observations = [
        observation_id for xtalform_site_id in assignments for observation_id in assignments[xtalform_site_id]
    ]
    observations_to_assign = [observation_id for observation_id in centroid_ca_positions]

    while len(observations_to_assign) > 0:
        for observation_id in observations_to_assign:
            # Assign observations near current xtalform sites
            for xtalform_site_id, xtalform_site_pos in centre_residues_positions.items():
                if _get_dist(centroid_ca_positions[observation_id], xtalform_site_pos) < cutoff:
                    assignments[xtalform_site_id].append(observation_id)
                    assigned_observations.append(observation_id)
                    break

        # Identify observation with most remaining items within 10A
        remaining_observations = [_k for _k in observations_to_assign if _k not in assigned_observations]
        if len(remaining_observations) == 0:
            break
        num_near = {
            observation_id: len(
                [
                    other_observation_id
                    for other_observation_id in remaining_observations
                    if (
                        _get_dist(centroid_ca_positions[observation_id], centroid_ca_positions[other_observation_id])
                        < cutoff
                    )
                    & (observation_id != other_observation_id)
                ]
            )
            for observation_id in remaining_observations
        }
        new_centroid_observation_id = min(num_near, key=lambda _x: num_near[_x])

        # Add it to centre residue list
        centre_residues_positions[new_centroid_observation_id] = centroid_ca_positions[new_centroid_observation_id]
        assignments[new_centroid_observation_id] = []

        # Update observations to assign
        # assigned_observations.append(new_centroid_observation_id)
        observations_to_assign = [_k for _k in observations_to_assign if _k not in assigned_observations]

    return assignments


def _update_xtalform_sites(
    xtalform_sites: dict[str, dt.XtalFormSite],
    canonical_site: dt.CanonicalSite,
    canonical_site_id: str,
    dataset_assignments: dict[str, str],
    conformer_sites: dict[str, dt.ConformerSite],
    neighbourhoods,
    debug=False,
):
    # Iterate over Canonical Sites, collecting Observations that are in the same crystalform
    # Then get their centroid CA positions,
    # spatially cluster on these with a reasonably broad cutoff, including old items
    # Then name these collections according to their centremost member

    # Partition by canonical site and xtalform
    crystalform_observations = {
        xtalform_name: [
            member
            for conformer_site_id in canonical_site.conformer_site_ids
            for member in conformer_sites[conformer_site_id].members
            if dataset_assignments[member[0]] == xtalform_name
        ]
        for xtalform_name in set(dataset_assignments.values())
    }

    # Get Observation centroid CA names and positions
    crystalform_observation_centroids = {
        xtalform_name: {
            member: _get_centroid_res(
                [_x for _x in set([(aid[0], aid[1]) for aid in neighbourhoods[member].atoms])], neighbourhoods[member]
            )
            for member in crystalform_observations[xtalform_name]
        }
        for xtalform_name in crystalform_observations
    }

    # Spatially cluster
    crystalform_observation_cluster_assignments = {
        xtalform_name: _crystalform_incremental_cluster(
            crystalform_observation_centroids[xtalform_name],
            {
                xid: xs
                for xid, xs in xtalform_sites.items()
                if (xs.xtalform_id == xtalform_name) & (xs.canonical_site_id == canonical_site_id)
            },
            neighbourhoods,
        )
        for xtalform_name in crystalform_observations
    }

    # if debug:
    #     raise Exception

    # Create the xtalforms or assign new observations
    for xtalform_name in crystalform_observation_cluster_assignments:
        for centroid_residue, asigned_observation_ids in crystalform_observation_cluster_assignments[
            xtalform_name
        ].items():
            # If the centroid is known, assign any new observations
            if centroid_residue in xtalform_sites:
                for asigned_observation_id in asigned_observation_ids:
                    if asigned_observation_id not in xtalform_sites[centroid_residue].members:
                        xtalform_sites[centroid_residue].members.append(asigned_observation_id)

            # Otherwise create a new crystalform site
            else:
                xtalform_site_id = "/".join(centroid_residue)
                xtalform_site = dt.XtalFormSite(
                    xtalform_name,
                    crystalform_observation_centroids[xtalform_name][centroid_residue][0],
                    canonical_site_id,
                    asigned_observation_ids,
                )
                xtalform_sites[xtalform_site_id] = xtalform_site

    # raise Exception
    #
    # matched = False
    # # for xtalform_site_id, xtalform_site in xtalform_sites.items():
    # #     if xtalform_site.canonical_site_id == canonical_site_id:
    # #         for conformer_site_id in canonical_site.conformer_site_ids:
    # #             conformer_site = conformer_sites[conformer_site_id]
    # #             for member in conformer_site.members:
    # #                 if dataset_assignments[member[0]] == xtalform_site.xtalform_id:
    # #                     if member not in xtalform_site.members:
    # #                         xtalform_site.members.append(member)
    #
    # xtalforms_dict = {
    #     (xtalform_site.canonical_site_id, xtalform_site.xtalform_id): xtalform_site_id
    #     for xtalform_site_id, xtalform_site in xtalform_sites.items()
    # }
    #
    # for conformer_site_id in canonical_site.conformer_site_ids:
    #     conformer_site = conformer_sites[conformer_site_id]
    #     for member in conformer_site.members:
    #         assignment = dataset_assignments[member[0]]
    #         if (canonical_site_id, assignment) in xtalforms_dict:
    #             xtalform_site = xtalform_sites[xtalforms_dict[(canonical_site_id, assignment)]]
    #             if member not in xtalform_site.members:
    #                 xtalform_site.members.append(member)
    #         else:
    #             xtalform_site_id = "/".join(member)
    #             xtalform_site = dt.XtalFormSite(
    #                 assignment,
    #                 member[1],
    #                 canonical_site_id,
    #                 [
    #                     member,
    #                 ],
    #             )
    #             xtalform_sites[xtalform_site_id] = xtalform_site
    #             xtalforms_dict[(xtalform_site.canonical_site_id, xtalform_site.xtalform_id)] = xtalform_site_id

    # Otherwise if not matched create a new xtalform site
    ...


def _save_xtalform_sites(fs_model, xtalform_sites: dict[str, dt.XtalFormSite]):
    with open(fs_model.xtalform_sites, "w") as f:
        dic = {}
        for xtalform_site_id, xtalform_site in xtalform_sites.items():
            dic[xtalform_site_id] = xtalform_site.to_dict()
        yaml.safe_dump(dic, f)


# def _update_conformer_site_transform(
#                 conformer_site_transforms,
#                 canonical_site,
#                 conformer_site,
#             ):
#
#
#     ...


def _save_conformer_site_transforms(
    fs_model: dt.FSModel, conformer_site_transforms: dict[tuple[str, str], dt.Transform]
):
    with open(fs_model.conformer_site_transforms, "w") as f:
        dic = {}
        for conformer_site_transform_id, conformer_site_transform in conformer_site_transforms.items():
            dic["~".join(conformer_site_transform_id)] = conformer_site_transform.to_dict()
        yaml.safe_dump(dic, f)
    ...


# def _update_canonical_site_transforms(
#             canonical_site_transforms: dict[tuple[str,str], dt.Transform],
#             canonical_site: dt.CanonicalSite,
#             canonical_sites: dict[str, dt.CanonicalSite],
#         ):

# ...


def _save_canonical_site_transforms(fs_model: dt.FSModel, canonical_site_transforms: dict[str, dt.Transform]):
    with open(fs_model.canonical_site_transforms, "w") as f:
        dic = {}
        for canonical_site_transform_id, canonical_site_transform in canonical_site_transforms.items():
            dic[canonical_site_transform_id] = canonical_site_transform.to_dict()
        yaml.safe_dump(dic, f)


def _update_fs_model(
    fs_model: dt.FSModel,
    canonical_sites: dict[str, dt.CanonicalSite],
    conformer_sites: dict[str, dt.ConformerSite],
    reference_datasets: dict[str, dt.Dataset],
    reference_transforms,
):
    # Iterate over canonical sites and their members, checking if they already have an output record and
    # if not creating one
    alignments = fs_model.alignments
    for canonical_site_id, canonical_site in canonical_sites.items():
        for conformer_site_id in canonical_site.conformer_site_ids:
            conformer_site = conformer_sites[conformer_site_id]
            for member in conformer_site.members:
                dtag, chain, residue, version = member
                if dtag not in alignments:
                    alignments[dtag] = {}
                if chain not in alignments[dtag]:
                    alignments[dtag][chain] = {}
                if residue not in alignments[dtag][chain]:
                    alignments[dtag][chain][residue] = {}
                if version not in alignments[dtag][chain][residue]:
                    alignments[dtag][chain][residue][version] = dt.LigandNeighbourhoodOutput(
                        {}, {}, {}, {}, {}, {}, {}, {}
                    )

                ligand_neighbourhood_output: dt.LigandNeighbourhoodOutput = alignments[dtag][chain][residue][version]

                if not (fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                    os.mkdir(fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag)

                if canonical_site_id not in ligand_neighbourhood_output.aligned_structures:
                    ligand_neighbourhood_output.aligned_structures[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_artefacts[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_xmaps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_XMAP_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_diff_maps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_DIFF_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_event_maps[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    # Crystallographic maps
                    ligand_neighbourhood_output.aligned_xmaps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_XMAP_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_diff_maps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_DIFF_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

                    ligand_neighbourhood_output.aligned_event_maps_crystallographic[canonical_site_id] = (
                        fs_model.source_dir
                        / constants.ALIGNED_FILES_DIR
                        / dtag
                        / constants.ALIGNED_EVENT_MAP_CRYSTALLOGRAPHIC_TEMPLATE.format(
                            dtag=dtag, chain=chain, residue=residue, version=version, site=canonical_site_id
                        )
                    )

    reference_alignments = fs_model.reference_alignments
    for dtag, dataset in reference_datasets.items():
        for canonical_site_id, canonical_site in canonical_sites.items():
            if (dtag, canonical_site_id) not in reference_transforms:
                continue

            if not (fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                os.mkdir(fs_model.source_dir / constants.ALIGNED_FILES_DIR / dtag)

            if dtag not in reference_alignments:
                reference_alignments[dtag] = {}

            if canonical_site_id not in reference_alignments[dtag]:
                reference_alignments[dtag][canonical_site_id] = {
                    "aligned_structures": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_STRUCTURE_TEMPLATE.format(dtag=dtag, site=canonical_site_id),
                    "aligned_artefacts": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                        dtag=dtag, site=canonical_site_id
                    ),
                    "aligned_xmaps": fs_model.source_dir
                    / constants.ALIGNED_FILES_DIR
                    / dtag
                    / constants.ALIGNED_REFERENCE_XMAP_TEMPLATE.format(dtag=dtag, site=canonical_site_id),
                    # 'aligned_event_maps': fs_model.source_dir / constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                    # dtag=dtag, chain=chain, residue=residue, site=canonical_site_id
                    # ),
                }


def _save_fs_model(fs_model: dt.FSModel):
    with open(fs_model.fs_model, "w") as f:
        dic = fs_model.to_dict()

        yaml.safe_dump(dic, f)


def save_reference_structure_transforms(
    fs_model: dt.FSModel,
    reference_structure_transforms: dict[tuple[str, str], dt.Transform],
):
    dic = {}
    for reference_structure_transform_id, reference_structure_transform in reference_structure_transforms.items():
        dic["~".join(reference_structure_transform_id)] = reference_structure_transform.to_dict()
    with open(fs_model.reference_structure_transforms, "w") as f:
        yaml.safe_dump(dic, f)


def _update(
    fs_model: dt.FSModel,
    datasets: dict[str, dt.Dataset],
    reference_datasets: dict[str, dt.Dataset],
    new_datasets: dict[str, dt.Dataset],
    assemblies: dict[str, dt.Assembly],
    xtalforms: dict[str, dt.XtalForm],
    dataset_assignments: dict[str, str],
    ligand_neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
    # alignment_landmarks: dict[tuple[str,str,str,int], dict[tuple[str, str, str], dt.Atom]],
    alignability_graph,
    connected_components,
    ligand_neighbourhood_transforms: dict[tuple[tuple[str, str, str, str], tuple[str, str, str, str]], dt.Transform],
    conformer_sites: dict[str, dt.ConformerSite],
    conformer_site_transforms: dict[tuple[str, str], dt.Transform],
    canonical_sites: dict[str, dt.CanonicalSite],
    # canonical_site_transforms: dict[str, dt.Transform],
    xtalform_sites: dict[str, dt.XtalFormSite],
    reference_structure_transforms: dict[tuple[str, str], dt.Transform],
    assembly_landmarks,
    assembly_transforms,
    version,
):
    logger.info(f"Version is: {version}")
    # Get the structures
    structures: dict = _get_structures(datasets)

    # Get the assembly alignment hierarchy
    hierarchy, biochain_priorities = alignment_heirarchy._derive_alignment_heirarchy(assemblies)
    alignment_heirarchy.save_yaml(fs_model.hierarchy, hierarchy, lambda x: x)
    alignment_heirarchy.save_yaml(fs_model.biochain_priorities, biochain_priorities, lambda x: x)

    # Get the assembly hierarchy transforms
    for assembly_name, assembly in assemblies.items():
        # Do not update if already have landmarks!
        if assembly_name in assembly_landmarks:
            continue
        as_st = alignment_heirarchy._get_assembly_st(assembly, structures[assembly.reference])
        assembly_landmarks[assembly_name] = alignment_heirarchy.structure_to_landmarks(as_st)
    alignment_heirarchy.save_yaml(
        fs_model.assembly_landmarks, assembly_landmarks, alignment_heirarchy.assembly_landmarks_to_dict
    )

    for assembly_name, assembly in assemblies.items():
        # Do not update if already have assembly transform!
        if assembly_name in assembly_transforms:
            continue
        assembly_transforms[assembly_name] = alignment_heirarchy._calculate_assembly_transform_sequence(
            hierarchy,
            assembly_name,
            assembly_landmarks,
        )
    alignment_heirarchy.save_yaml(fs_model.assembly_transforms, assembly_transforms, lambda x: x)

    # Assign datasets
    new_dataset_assignments = {}
    for dtag, dataset in new_datasets.items():
        new_dataset_assignments[dtag] = _assign_dataset(
            dataset,
            assemblies,
            xtalforms,
            structures[dtag],
            structures,
        )
    dataset_assignments.update(new_dataset_assignments)

    _save_assignments(fs_model, dataset_assignments)
    logger.info(f"Assigned {len(dataset_assignments)} xtalform assignments to datasets!")

    # Get neighbourhoods
    logger.info(f"Updating neighbourhoods")
    for dtag, dataset in new_datasets.items():
        xtalform = xtalforms[dataset_assignments[dtag]]
        neighborhoods = _get_neighbourhoods(dataset, xtalform, assemblies, version)
        logger.info(f"Dataset {dtag} has {len(neighborhoods)} ligand neighbourhoods")
        for lid, neighbourhood in neighborhoods.items():
            ligand_neighbourhoods[lid] = neighbourhood
    logger.info(f"Found {len(ligand_neighbourhoods)} ligand neighbourhoods!")
    _save_neighbourhoods(fs_model, ligand_neighbourhoods)

    # for nid, neighbourhood in ligand_neighbourhoods.items():
    #     print(nid)
    #     for atom_id in neighbourhood.atoms:
    #         if atom_id[2] == "CA":
    #             print(atom_id)

    # Get chain to assembly transforms
    logger.info(f"Getting chain-to-assembly transforms...")
    chain_to_assembly_transforms = {}
    for dtag, st in structures.items():
        xtalform_chains = [
            _chain
            for _xassembly in xtalforms[dataset_assignments[dtag]].assemblies.values()
            for _chain in _xassembly.chains
        ]

        dataset_chains = [_chain.name for _chain in st[0]]
        dataset_ligand_chains = [_x[1] for _x in ligand_neighbourhoods if _x[0] == dtag]
        # for _chain in dataset_ligand_chains:
        # for _chain in dataset_chains:
        for _chain in xtalform_chains:
            if _chain not in xtalform_chains:
                raise Exception(
                    f"A xtalform assignment error has occured. Dataset {dtag} has chain {_chain} in its chains {dataset_chains} however its assigned xtalform {dataset_assignments[dtag]} has chain {xtalform_chains}"
                )
            try:
                chain_to_assembly_transforms[
                    (
                        dtag,
                        _chain,
                        # version,
                    )
                ] = alignment_heirarchy._get_structure_chain_to_assembly_transform(
                    st,
                    _chain,
                    xtalforms[dataset_assignments[dtag]],
                    assemblies,
                    assembly_landmarks,
                )
            except Exception as e:
                print(f'Exception in dataset: {dtag}, in xtalform {dataset_assignments[dtag]}')
                raise e
    logger.info(f'Got {len(chain_to_assembly_transforms)} chain to assembly transforms')
    alignment_heirarchy.save_yaml(
        fs_model.chain_to_assembly,
        chain_to_assembly_transforms,
        alignment_heirarchy.chain_to_assembly_transforms_to_dict,
    )
    logger.info(f'Got {len(chain_to_assembly_transforms)} chain to assembly transforms')

    # Update graph
    logger.info(f"Updating alignment graph...")
    logger.info(f"Previously had {len(ligand_neighbourhood_transforms)} alignments between neighbourhoods")

    # for dtag, dataset in new_datasets.items():
    for lid, neighbourhood in ligand_neighbourhoods.items():
        _update_ligand_neighbourhood_transforms(
            ligand_neighbourhood_transforms,
            lid,
            ligand_neighbourhoods,
            structures,
        )

    logger.info(f"Now have {len(ligand_neighbourhood_transforms)} alignments between neighbourhoods")
    # print(ligand_neighbourhood_transforms)
    _save_ligand_neighbourhood_transforms(fs_model, ligand_neighbourhood_transforms)

    # Update the alignment graph
    logger.info(f"Updating alignment graph...")
    logger.info(f"Previously had {len(alignability_graph.nodes)} nodes")
    logger.info(f"Previously had {len(alignability_graph.edges)} edges")
    _update_graph(
        alignability_graph,
        ligand_neighbourhoods,
        ligand_neighbourhood_transforms,
    )
    logger.info(f"Now have {len(alignability_graph.nodes)} nodes")
    # print(sorted(alignability_graph.nodes, key=lambda x: x[0]))
    logger.info(f"Now have {len(alignability_graph.edges)} edges")
    # print(alignability_graph.edges)
    _save_graph(fs_model, alignability_graph)

    # Update conformer sites
    logger.info(f"Updating conformer sites...")
    connected_components = _get_connected_components(alignability_graph, connected_components)
    _save_connected_components(fs_model, connected_components)
    logger.info(f"Got {len(connected_components)} connected components")
    logger.info(f"Previously had {len(conformer_sites)} conformer sites")

    for connected_component_id, connected_component in connected_components.items():
        # Update new datasets to indicate everything sharing a connected component
        #

        # Match new component to old ones by membership, and expand old ones if available otherwise create new one
        _update_conformer_sites(
            conformer_sites,
            connected_component_id,
            connected_component,
            ligand_neighbourhoods,
            structures,
            assemblies,
            xtalforms,
            dataset_assignments,
        )
    logger.info(f"Now have {len(conformer_sites)} conformer sites")
    # for conformer_site_id, conformer_site in conformer_sites.items():
    #     print(conformer_site_id)
    #     print(conformer_site.members)
    _save_conformer_sites(fs_model, conformer_sites)

    # Update canonical sites
    logger.info(f"Previously had {len(canonical_sites)} canonical sites")
    for conformer_site_id, conformer_site in conformer_sites.items():
        # If conformer site in a canonical site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        _update_canonical_sites(canonical_sites, conformer_site, conformer_site_id, ligand_neighbourhoods)
    logger.info(f"Now have {len(canonical_sites)} canonical sites")
    logger.info(f"Global reference dtag is: {list(canonical_sites.values())[0].global_reference_dtag}")
    _save_canonical_sites(fs_model, canonical_sites)

    # Update crystalform sites
    logger.info(f"Previously had {len(xtalform_sites)} xtalform sites")
    for canonical_site_id, canonical_site in canonical_sites.items():
        # If canonical site in a xtalform site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        debug = False
        if canonical_site_id == 'LYSRSCPZ-x0426+A+805+1':
            debug = True
        _update_xtalform_sites(
            xtalform_sites,
            canonical_site,
            canonical_site_id,
            dataset_assignments,
            conformer_sites,
            ligand_neighbourhoods,
            debug,
        )
    logger.info(f"Now have {len(xtalform_sites)} xtalform sites")
    _save_xtalform_sites(fs_model, xtalform_sites)

    # Get conformer site transforms
    logger.info(f"Previously had {len(conformer_site_transforms)} conformer site transforms")
    for canonical_site_id, canonical_site in canonical_sites.items():
        # for conformer_site_id in canonical_site.conformer_site_ids:
        # conformer_site = conformer_sites[conformer_site_id]
        # print(conformer_site)
        _update_conformer_site_transforms(
            conformer_site_transforms,
            canonical_site,
            conformer_sites,
            structures,
        )
    logger.info(f"Now have {len(conformer_site_transforms)} conformer site transforms")
    # for conformer_site_transform_id, conformer_site_transform in conformer_site_transforms.items():
    #     print(conformer_site_transform_id)
    _save_conformer_site_transforms(fs_model, conformer_site_transforms)

    # Get canonical site tranforms
    # logger.info(f"Previously had {len(canonical_site_transforms)} canonical site transforms...")
    # for canonical_site_id, canonical_site in canonical_sites.items():
    #     _update_canonical_site_transforms(
    #         canonical_site_transforms,
    #         canonical_site_id,
    #         canonical_site,
    #         conformer_sites,
    #         structures
    #     )
    #     logger.info(f"Now have {len(canonical_site_transforms)} canonical site transforms")
    # _save_canonical_site_transforms(fs_model, canonical_site_transforms)

    # Update the reference structure transforms
    logger.info(f"Previously had {len(reference_structure_transforms)} reference structure transforms")
    for dtag, dataset in reference_datasets.items():
        for canonical_site_id, canonical_site in canonical_sites.items():
            key = (dtag, canonical_site_id)
            if key not in reference_structure_transforms:
                _update_reference_structure_transforms(
                    reference_structure_transforms,
                    key,
                    structures,
                    canonical_site,
                    conformer_sites,
                    assemblies,
                    xtalforms,
                    dataset_assignments,
                    xtalform_sites,
                    canonical_site_id,
                )
    logger.info(f"Now have {len(reference_structure_transforms)} reference structure transforms")
    save_reference_structure_transforms(
        fs_model,
        reference_structure_transforms,
    )

    # Update output: check if aligned data for each lid in canon site is already there and if not add it
    _update_fs_model(fs_model, canonical_sites, conformer_sites, reference_datasets, reference_structure_transforms)
    _save_fs_model(fs_model)

    # Generate new aligned structures
    # for canonical_site_id, canonical_site in canonical_sites.items():
    #     for conformer_site_id, conformer_site in canonical_site.conformer_sites.items():
    #         for lid in conformer_site.ligand_ids:

    for dtag, dataset_alignment_info in fs_model.alignments.items():
        for chain, chain_alignment_info in dataset_alignment_info.items():
            for residue, residue_alignment_info in chain_alignment_info.items():
                for version, ligand_neighbourhood_output in residue_alignment_info.items():
                    for (
                        canonical_site_id,
                        aligned_structure_path,
                    ) in ligand_neighbourhood_output.aligned_structures.items():
                        if not (
                            (fs_model.source_dir.parent / aligned_structure_path).exists()
                            | Path(aligned_structure_path).exists()
                        ):
                            # _update_aligned_structures()
                            _structure = structures[dtag].clone()
                            canonical_site = canonical_sites[canonical_site_id]
                            # Check for the matching conformer site
                            conformer_site = None
                            for conformer_site_id in canonical_site.conformer_site_ids:
                                if (dtag, chain, residue, version) in conformer_sites[conformer_site_id].members:
                                    conformer_site = conformer_sites[conformer_site_id]
                                    break
                            if conformer_site is None:
                                print(f"Skipping alignment of {dtag} {chain} {residue} to site {canonical_site_id}!")
                                continue
                            moving_ligand_id = (dtag, chain, residue, version)
                            reference_ligand_id = conformer_site.reference_ligand_id
                            # print(aligned_structure_path)

                            # Get the site chain
                            site_chain = None
                            xtalform_site = None
                            # site_reference_ligand_id = conformer_sites[
                            #     canonical_site.reference_conformer_site_id].reference_ligand_id
                            site_reference_ligand_id = conformer_site.reference_ligand_id
                            site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
                            site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]
                            for xsid, _xtalform_site in xtalform_sites.items():
                                _xtalform_id = _xtalform_site.xtalform_id
                                _xtalform_canonical_site_id = _xtalform_site.canonical_site_id
                                if (
                                    (_xtalform_id == site_reference_ligand_xtalform_id)
                                    & (_xtalform_canonical_site_id == canonical_site_id)
                                    & (site_reference_ligand_id in _xtalform_site.members)
                                ):
                                    xtalform_site = _xtalform_site
                            site_chain = xtalform_site.crystallographic_chain

                            # Aligns to conformer site, then to the corresponding assembly, then from that assembly to
                            # the global frame
                            try:
                                _align_structure(
                                    _structure,
                                    moving_ligand_id,
                                    reference_ligand_id,
                                    ligand_neighbourhoods[moving_ligand_id],
                                    [nid for nid in ligand_neighbourhoods if nid[0] == dtag],
                                    alignability_graph,
                                    ligand_neighbourhood_transforms,
                                    conformer_site_transforms,
                                    # canonical_site_transforms,
                                    canonical_site_id,
                                    conformer_site_id,
                                    xtalforms[dataset_assignments[dtag]],
                                    aligned_structure_path,
                                    site_reference_xform=xtalforms[
                                        dataset_assignments[conformer_site.reference_ligand_id[0]]
                                    ],
                                    chain_to_assembly_transform=chain_to_assembly_transforms[
                                        (
                                            conformer_site.reference_ligand_id[0],
                                            # conformer_site.reference_ligand_id[1],
                                            site_chain
                                            # conformer_site.reference_ligand_id[3]
                                        )
                                    ],
                                    assembly_transform=assembly_transforms[
                                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                        .assemblies[
                                            alignment_heirarchy._chain_to_xtalform_assembly(
                                                # conformer_site.reference_ligand_id[1],
                                                site_chain,
                                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                                            )
                                        ]
                                        .assembly
                                    ],
                                    xtalform_sites=xtalform_sites,
                                )

                            except:
                                logger.info(f"Failed to generate aligned structure {aligned_structure_path}")
                                raise
                        else:
                            logger.info(f"Already output structure!")
    # Generate alignments of references to each canonical site
    # for canonical_site_id, canonical_site in canonical_sites.items():
    #     for dtag, reference_dataset in reference_datasets.items():
    # _update_reference_alignments(
    #
    # )
    for dtag, dataset_alignment_info in fs_model.reference_alignments.items():
        for canonical_site_id, alignment_info in dataset_alignment_info.items():
            aligned_structure_path = alignment_info["aligned_structures"]
            logger.info(f"Outputting reference structure: {aligned_structure_path}")
            if not (
                (fs_model.source_dir.parent / aligned_structure_path).exists() | Path(aligned_structure_path).exists()
            ):
                _structure = structures[dtag].clone()
                _align_reference_structure(
                    _structure,
                    dtag,
                    reference_structure_transforms,
                    # canonical_site_transforms,
                    canonical_site_id,
                    alignment_info["aligned_structures"],
                )
            else:
                logger.info(f"Already output reference structure!")

    # Generate new aligned maps
    # for canonical_site_id, canonical_site in canonical_sites.items():
    #     for conformer_site_id in canonical_site.conformer_site_ids:
    #         conformer_site = conformer_sites[conformer_site_id]
    #         for lid in conformer_site.members:
    #             _update_aligned_xmaps()
    reference_xmap = read_xmap_from_mtz(datasets[[x for x in canonical_sites.values()][0].global_reference_dtag].mtz)
    logger.info(f"Outputting xmaps...")
    for dtag, dataset_alignment_info in fs_model.alignments.items():
        for chain, chain_alignment_info in dataset_alignment_info.items():
            for residue, residue_alignment_info in chain_alignment_info.items():
                for version, ligand_neighbourhood_output in residue_alignment_info.items():
                    for (
                        canonical_site_id,
                        aligned_event_map_path,
                    ) in ligand_neighbourhood_output.aligned_event_maps.items():
                        logger.info(f"Writing to: {aligned_event_map_path}")
                        # if not ( (fs_model.source_dir.parent / aligned_event_map_path).exists() | Path(aligned_event_map_path).exists()):
                        rel = True
                        try:
                            rel = Path(aligned_event_map_path).is_relative_to(fs_model.source_dir)
                        except:
                            rel = False
                        if rel:
                            # print(fs_model.source_dir.parent / aligned_event_map_path)
                            _structure = structures[dtag].clone()
                            canonical_site = canonical_sites[canonical_site_id]
                            # Check for the matching conformer site
                            conformer_site = None
                            for conformer_site_id in canonical_site.conformer_site_ids:
                                if (dtag, chain, residue, version) in conformer_sites[conformer_site_id].members:
                                    conformer_site = conformer_sites[conformer_site_id]
                                    break

                            if conformer_site is None:
                                print(f"Skipping alignment of {dtag} {chain} {residue} to site {canonical_site_id}!")
                                continue

                            moving_ligand_id = (dtag, chain, residue, version)
                            reference_ligand_id = conformer_site.reference_ligand_id
                            # print(ligand_neighbourhoods)

                            xmap_path = datasets[dtag].ligand_binding_events[(dtag, chain, residue)].xmap

                            aligned_structure_path = ligand_neighbourhood_output.aligned_structures[canonical_site_id]
                            st_path = fs_model.source_dir.parent / aligned_structure_path
                            if not st_path.exists():
                                st_path = aligned_structure_path
                            aligned_structure = gemmi.read_structure(str(st_path))
                            aligned_res = aligned_structure[0][chain][str(residue)][0]

                            # Get the site chain
                            site_chain = None
                            xtalform_site = None
                            # site_reference_ligand_id = conformer_sites[
                            #     canonical_site.reference_conformer_site_id].reference_ligand_id
                            site_reference_ligand_id = conformer_site.reference_ligand_id
                            site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
                            site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]
                            for xsid, _xtalform_site in xtalform_sites.items():
                                _xtalform_id = _xtalform_site.xtalform_id
                                _xtalform_canonical_site_id = _xtalform_site.canonical_site_id
                                if (
                                    (_xtalform_id == site_reference_ligand_xtalform_id)
                                    & (_xtalform_canonical_site_id == canonical_site_id)
                                    & (site_reference_ligand_id in _xtalform_site.members)
                                ):
                                    xtalform_site = _xtalform_site
                            site_chain = xtalform_site.crystallographic_chain

                            if (xmap_path != "None") and (xmap_path is not None):
                                xmap = read_xmap(xmap_path)

                                __align_xmap(
                                    ligand_neighbourhoods[(dtag, chain, residue, version)],
                                    alignability_graph,
                                    ligand_neighbourhood_transforms,
                                    reference_xmap,
                                    reference_ligand_id,
                                    moving_ligand_id,
                                    xmap,
                                    conformer_site_transforms,
                                    conformer_site_id,
                                    # canonical_site_transforms,
                                    canonical_site_id,
                                    aligned_event_map_path,
                                    ligand_neighbourhood_output.aligned_event_maps_crystallographic[canonical_site_id],
                                    aligned_res,
                                    chain_to_assembly_transform=chain_to_assembly_transforms[
                                        (
                                            conformer_site.reference_ligand_id[0],
                                            # conformer_site.reference_ligand_id[1],
                                            # alignment_heirarchy._chain_to_xtalform_assembly(
                                            #     # conformer_site.reference_ligand_id[1],
                                            #     site_chain,
                                            #     xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                            # )
                                            site_chain
                                            # conformer_site.reference_ligand_id[3]
                                        )
                                    ],
                                    assembly_transform=assembly_transforms[
                                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                        .assemblies[
                                            alignment_heirarchy._chain_to_xtalform_assembly(
                                                site_chain,
                                                # conformer_site.reference_ligand_id[1],
                                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                                            )
                                        ]
                                        .assembly
                                    ],
                                )
                            mtz_path = datasets[dtag].mtz
                            # print(f"Mtz path: {mtz_path}")
                            # raise Exception
                            if mtz_path != "None":
                                xmap = read_xmap_from_mtz(mtz_path, "2Fo-Fc")
                                __align_xmap(
                                    ligand_neighbourhoods[(dtag, chain, residue, version)],
                                    alignability_graph,
                                    ligand_neighbourhood_transforms,
                                    reference_xmap,
                                    reference_ligand_id,
                                    moving_ligand_id,
                                    xmap,
                                    conformer_site_transforms,
                                    conformer_site_id,
                                    # canonical_site_transforms,
                                    canonical_site_id,
                                    ligand_neighbourhood_output.aligned_xmaps[canonical_site_id],
                                    ligand_neighbourhood_output.aligned_xmaps_crystallographic[canonical_site_id],
                                    aligned_res,
                                    chain_to_assembly_transform=chain_to_assembly_transforms[
                                        (
                                            conformer_site.reference_ligand_id[0],
                                            # conformer_site.reference_ligand_id[1],
                                            site_chain
                                            # alignment_heirarchy._chain_to_xtalform_assembly(
                                            #     site_chain,
                                            #     # conformer_site.reference_ligand_id[1],
                                            #     xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                            # )
                                            # conformer_site.reference_ligand_id[3]
                                        )
                                    ],
                                    assembly_transform=assembly_transforms[
                                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                        .assemblies[
                                            alignment_heirarchy._chain_to_xtalform_assembly(
                                                # conformer_site.reference_ligand_id[1],
                                                site_chain,
                                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                                            )
                                        ]
                                        .assembly
                                    ],
                                )
                                xmap = read_xmap_from_mtz(mtz_path, "Fo-Fc")
                                __align_xmap(
                                    ligand_neighbourhoods[(dtag, chain, residue, version)],
                                    alignability_graph,
                                    ligand_neighbourhood_transforms,
                                    reference_xmap,
                                    reference_ligand_id,
                                    moving_ligand_id,
                                    xmap,
                                    conformer_site_transforms,
                                    conformer_site_id,
                                    # canonical_site_transforms,
                                    canonical_site_id,
                                    ligand_neighbourhood_output.aligned_diff_maps[canonical_site_id],
                                    ligand_neighbourhood_output.aligned_diff_maps_crystallographic[canonical_site_id],
                                    aligned_res,
                                    chain_to_assembly_transform=chain_to_assembly_transforms[
                                        (
                                            conformer_site.reference_ligand_id[0],
                                            # conformer_site.reference_ligand_id[1],
                                            site_chain
                                            # conformer_site.reference_ligand_id[3]
                                        )
                                    ],
                                    assembly_transform=assembly_transforms[
                                        xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]]
                                        .assemblies[
                                            alignment_heirarchy._chain_to_xtalform_assembly(
                                                site_chain,
                                                # conformer_site.reference_ligand_id[1],
                                                xtalforms[dataset_assignments[conformer_site.reference_ligand_id[0]]],
                                            )
                                        ]
                                        .assembly
                                    ],
                                )

                        else:
                            logger.info(f"Already output xmap!")
    return fs_model


def _load_assemblies(assemblies_file, new_assemblies_yaml):
    assemblies = {}

    if assemblies_file.exists():
        with open(assemblies_file, "r") as f:
            dic = yaml.safe_load(f)

        for assembly_id, assembly_info in dic["assemblies"].items():
            assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    # Load new info and update
    if new_assemblies_yaml.exists():
        with open(new_assemblies_yaml, "r") as f:
            new_assemblies_dict = yaml.safe_load(f)["assemblies"]
    else:
        new_assemblies_dict = {}

    for assembly_id, assembly_info in new_assemblies_dict.items():
        if assembly_id in assemblies:
            continue
        assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return assemblies


def _load_xtalforms(xtalforms_file, new_xtalforms_yaml):
    xtalforms = {}

    if xtalforms_file.exists():
        with open(xtalforms_file, "r") as f:
            dic = yaml.safe_load(f)["crystalforms"]

        for xtalform_id, xtalform_info in dic.items():
            xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    # Load new info and update
    if new_xtalforms_yaml.exists():
        with open(new_xtalforms_yaml, "r") as f:
            new_xtalforms_dict = yaml.safe_load(f)["crystalforms"]
    else:
        new_xtalforms_dict = {}

    for xtalform_id, xtalform_info in new_xtalforms_dict.items():
        if xtalform_id in xtalforms:
            continue
        xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return xtalforms


def _load_xtalforms_and_assemblies(xtalforms_file, new_xtalforms_yaml):
    assemblies = {}
    xtalforms = {}

    if xtalforms_file.exists():
        with open(xtalforms_file, "r") as f:
            dic = yaml.safe_load(f)

        for xtalform_id, xtalform_info in dic.items():
            xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    # Load new info and update
    if new_xtalforms_yaml.exists():
        with open(new_xtalforms_yaml, "r") as f:
            new_xtalforms_dict = yaml.safe_load(f)
    else:
        new_xtalforms_dict = {}

    for xtalform_id, xtalform_info in new_xtalforms_dict.items():
        if xtalform_id in xtalforms:
            continue
        xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return xtalforms


def _load_dataset_assignments(dataset_assignments_yaml, fail_if_not_found=False):
    if fail_if_not_found and not dataset_assignments_yaml.is_file():
        raise ValueError("File " + str(dataset_assignments_yaml) + " does not exist")
    dataset_assignments = {}
    if dataset_assignments_yaml.exists():
        with open(dataset_assignments_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for dtag, assignment in dic.items():
            dataset_assignments[dtag] = assignment

    return dataset_assignments


def _load_ligand_neighbourhoods(ligand_neighbourhoods_yaml, fail_if_not_found=False):
    if fail_if_not_found and not ligand_neighbourhoods_yaml.is_file():
        raise ValueError("File " + str(ligand_neighbourhoods_yaml) + " does not exist")
    ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = {}

    if ligand_neighbourhoods_yaml.exists():
        with open(ligand_neighbourhoods_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic:
            for ligand_id, neighbourhood_info in dic.items():
                dtag, chain, residue, version = ligand_id.split("/")
                neighbourhood = dt.Neighbourhood.from_dict(neighbourhood_info)
                ligand_neighbourhoods[(dtag, chain, residue, version)] = neighbourhood

    return ligand_neighbourhoods


def _load_alignability_graph(alignability_graph, fail_if_not_found=False):
    if fail_if_not_found and not alignability_graph.is_file():
        raise ValueError("File " + str(alignability_graph) + " does not exist")
    if alignability_graph.exists():
        g_initial = nx.read_gml(
            str(alignability_graph),
            # destringizer=lambda x: tuple(x.split("/")),
        )

        g = nx.relabel_nodes(g_initial, {x: tuple(x.split("/")) for x in g_initial})

        return g

    else:
        return nx.Graph()


def _load_connected_components(connected_components_yaml, fail_if_not_found=False):
    if fail_if_not_found and not connected_components_yaml.is_file():
        raise ValueError("File " + str(connected_components_yaml) + " does not exist")
    connected_components = {}

    if connected_components_yaml.exists():
        with open(connected_components_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic:
            for ligand_id, neighbourhood_info in dic.items():
                dtag, chain, residue, version = ligand_id.split("+")
                # neighbourhood = dt.Neighbourhood.from_dict(neighbourhood_info)
                connected_components[(dtag, chain, residue, version)] = [
                    tuple([x for x in _ligand_id.split("+")]) for _ligand_id in neighbourhood_info
                ]

    return connected_components


def _load_ligand_neighbourhood_transforms(ligand_neighbourhood_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not ligand_neighbourhood_transforms_yaml.is_file():
        raise ValueError("File " + str(ligand_neighbourhood_transforms_yaml) + " does not exist")
    ligand_neighbourhood_transforms = {}
    if ligand_neighbourhood_transforms_yaml.exists():
        with open(ligand_neighbourhood_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for ligand_transform_key, ligand_transform in dic.items():
            # print(ligand_transform_key)
            ligand_1_id, ligand_2_id = ligand_transform_key.split("~")
            dtag_1, chain_1, residue_1, version_1 = ligand_1_id.split("/")
            dtag_2, chain_2, residue_2, version_2 = ligand_2_id.split("/")
            ligand_neighbourhood_transforms[
                ((dtag_1, chain_1, residue_1, version_1), (dtag_2, chain_2, residue_2, version_2))
            ] = dt.Transform.from_dict(ligand_transform)

    return ligand_neighbourhood_transforms


def _load_conformer_sites(conformer_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not conformer_sites_yaml.is_file():
        raise ValueError("File " + str(conformer_sites_yaml) + " does not exist")
    conformer_sites = {}
    if conformer_sites_yaml.exists():
        with open(conformer_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)
        for conformer_site_id, conformer_site_info in dic.items():
            conformer_sites[conformer_site_id] = dt.ConformerSite.from_dict(conformer_site_info)

    return conformer_sites


def _load_conformer_site_transforms(conformer_site_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not conformer_site_transforms_yaml.is_file():
        raise ValueError("File " + str(conformer_site_transforms_yaml) + " does not exist")
    conformer_site_transforms = {}
    if conformer_site_transforms_yaml.exists():
        with open(conformer_site_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for conformer_site_transform_id, conformer_site_transform_info in dic.items():
            conformer_site_1, conformer_site_2 = conformer_site_transform_id.split("~")

            conformer_site_transforms[(conformer_site_1, conformer_site_2)] = dt.Transform.from_dict(
                conformer_site_transform_info
            )

    return conformer_site_transforms


def _load_canonical_sites(canonical_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not canonical_sites_yaml.is_file():
        raise ValueError("File " + str(canonical_sites_yaml) + " does not exist")
    canonical_sites = {}
    if canonical_sites_yaml.exists():
        with open(canonical_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic is not None:
            for canonical_site_id, canonical_site_info in dic.items():
                canonical_sites[canonical_site_id] = dt.CanonicalSite.from_dict(canonical_site_info)

    return canonical_sites


def _load_canonical_site_transforms(canonical_site_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not canonical_site_transforms_yaml.is_file():
        raise ValueError("File " + str(canonical_site_transforms_yaml) + " does not exist")
    canonical_site_transforms = {}
    if canonical_site_transforms_yaml.exists():
        with open(canonical_site_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for canonical_site_id, canonical_site_transform_info in dic.items():
            canonical_site_transforms[canonical_site_id] = dt.Transform.from_dict(canonical_site_transform_info)

    return canonical_site_transforms


def _load_xtalform_sites(xtalform_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not xtalform_sites_yaml.is_file():
        raise ValueError("File " + str(xtalform_sites_yaml) + " does not exist")
    xtalform_sites = {}
    if xtalform_sites_yaml.exists():
        with open(xtalform_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for xtalform_site_id, xtalform_site_info in dic.items():
            xtalform_sites[xtalform_site_id] = dt.XtalFormSite.from_dict(xtalform_site_info)

    return xtalform_sites


def _load_reference_stucture_transforms(reference_structure_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not reference_structure_transforms_yaml.is_file():
        raise ValueError("File " + str(reference_structure_transforms_yaml) + " does not exist")
    reference_structure_transforms = {}
    if reference_structure_transforms_yaml.exists():
        with open(reference_structure_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for reference_structure_transform_id, reference_structure_transform_info in dic.items():
            dtag, canonical_site_id = reference_structure_transform_id.split("~")
            reference_structure_transforms[(dtag, canonical_site_id)] = dt.Transform.from_dict(
                reference_structure_transform_info
            )

    return reference_structure_transforms


class CLI:
    def schema(self, output_dir: str):
        _output_dir = Path(output_dir)

        if not _output_dir.exists():
            os.mkdir(_output_dir)

        save_schema(SystemData, _output_dir)
        save_schema(LigandNeighbourhoods, _output_dir)
        save_schema(CanonicalSites, _output_dir)
        save_schema(Transforms, _output_dir)
        save_schema(SiteTransforms, _output_dir)
        save_schema(SystemData, _output_dir)
        save_schema(SystemData, _output_dir)
        save_schema(SystemData, _output_dir)

    def update(self, options_json: str):
        options = Options.parse_file(options_json)
        # logger.info(f"Options json path is: {options}")
        logger.info(f"Input dir is: {options.source_dir}")
        logger.info(f"Output dir is: {options.output_dir}")

        if options.source_dir:
            source_fs_model = dt.FSModel.from_dir(options.source_dir)
        else:
            source_fs_model = None
        # else:
        # source_fs_model = dt.FSModel.default()

        fs_model = dt.FSModel.from_dir(
            options.output_dir,
        )
        if source_fs_model:
            fs_model.alignments = source_fs_model.alignments
            fs_model.reference_alignments = source_fs_model.reference_alignments

        # Update the output fs model, creating flat symlinks to old data

        if not Path(options.output_dir).exists():
            os.mkdir(options.output_dir)

        aligned_files_dir = Path(options.output_dir) / constants.ALIGNED_FILES_DIR
        if not aligned_files_dir.exists():
            os.mkdir(aligned_files_dir)

        # aligned_xmap_dir = Path(options.output_dir) / constants.ALIGNED_XMAPS_DIR
        # if not aligned_xmap_dir.exists():
        #     os.mkdir(aligned_xmap_dir)
        fs_model.symlink_old_data()

        source_data_model = dt.SourceDataModel.from_fs_model(
            fs_model, options.datasources, options.datasource_types, options.panddas
        )

        datasets, reference_datasets, new_datasets = source_data_model.get_datasets()
        # print(datasets)
        # print(reference_datasets)
        # print(new_datasets)

        # Get assemblies
        logger.info(f"Getting assemblies...")
        if source_fs_model:
            assemblies: dict[str, dt.Assembly] = _load_assemblies(
                source_fs_model.assemblies, Path(options.assemblies_json)
            )
        else:
            assemblies = _load_assemblies(fs_model.assemblies, Path(options.assemblies_json))
        # for key, assembly in assemblies.items():
        #     print(assembly)
        #     for gen in assembly.generators:
        #         print([gen.chain, gen.reference_chain, gen.triplet])
        # print(assemblies)

        # Get xtalforms
        logger.info(f"Getting xtalforms...")
        if source_fs_model:
            xtalforms: dict[str, dt.XtalForm] = _load_xtalforms(
                source_fs_model.xtalforms, Path(options.xtalforms_json)
            )
        else:
            xtalforms = _load_xtalforms(fs_model.xtalforms, Path(options.xtalforms_json))
        # for key, xtalform in xtalforms.items():
        #     print(xtalform)
        #     for ass, xtalform_ass in xtalform.assemblies.items():
        #         for chn, trns in zip(xtalform_ass.chains, xtalform_ass.transforms):
        #             print([chn, trns])
        # print(xtalforms)

        # Get the dataset assignments
        logger.info(f"Getting dataset assignments...")
        if source_fs_model:
            dataset_assignments = _load_dataset_assignments(Path(source_fs_model.dataset_assignments))
        else:
            dataset_assignments = _load_dataset_assignments(Path(fs_model.dataset_assignments))
        # print(dataset_assignments)

        # Get Ligand neighbourhoods
        logger.info(f"Getting ligand neighbourhoods...")
        if source_fs_model:
            ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = _load_ligand_neighbourhoods(
                source_fs_model.ligand_neighbourhoods
            )
        else:
            ligand_neighbourhoods = _load_ligand_neighbourhoods(fs_model.ligand_neighbourhoods)
        # print(ligand_neighbourhoods)

        # Get alignability graph
        logger.info(f"Getting alignability graph...")
        if source_fs_model:
            alignability_graph = _load_alignability_graph(source_fs_model.alignability_graph)
        else:
            alignability_graph = _load_alignability_graph(fs_model.alignability_graph)

        if source_fs_model:
            connected_components = _load_connected_components(source_fs_model.connected_components)
        else:
            connected_components = _load_connected_components(fs_model.connected_components)

        #
        logger.info(f"Getting lighand neighbourhood transforms...")
        if source_fs_model:
            ligand_neighbourhood_transforms: dict[
                tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform
            ] = _load_ligand_neighbourhood_transforms(source_fs_model.ligand_neighbourhood_transforms)
        else:
            ligand_neighbourhood_transforms = _load_ligand_neighbourhood_transforms(
                fs_model.ligand_neighbourhood_transforms
            )

        # Get conformer sites
        logger.info(f"Getting conformer sites...")
        if source_fs_model:
            conformer_sites: dict[str, dt.ConformerSite] = _load_conformer_sites(source_fs_model.conformer_sites)
        else:
            conformer_sites = _load_conformer_sites(fs_model.conformer_sites)

        #
        logger.info(f"Getting conformer site transforms...")
        if source_fs_model:
            conformer_site_transforms: dict[tuple[str, str], dt.Transform] = _load_conformer_site_transforms(
                source_fs_model.conformer_site_transforms
            )
        else:
            conformer_site_transforms = _load_conformer_site_transforms(fs_model.conformer_site_transforms)

        # Get canonical sites
        logger.info(f"Getting canonical sites...")
        if source_fs_model:
            canonical_sites: dict[str, dt.CanonicalSite] = _load_canonical_sites(source_fs_model.canonical_sites)
        else:
            canonical_sites = _load_canonical_sites(fs_model.canonical_sites)

        #
        # logger.info(f"Getting canonical site transforms...")
        # if source_fs_model:
        #     canonical_site_transforms: dict[str, dt.Transform] = _load_canonical_site_transforms(
        #         source_fs_model.conformer_site_transforms)
        # else:
        #     canonical_site_transforms = _load_canonical_site_transforms(
        #         fs_model.conformer_site_transforms)

        # Get xtalform sites
        logger.info(f"Getting xtalform sites...")
        if source_fs_model:
            xtalform_sites: dict[str, dt.XtalFormSite] = _load_xtalform_sites(source_fs_model.xtalform_sites)
        else:
            xtalform_sites = _load_xtalform_sites(fs_model.xtalform_sites)

        # Get reference structure transforms
        logger.info(f"Getting reference structure transforms...")
        if source_fs_model:
            reference_structure_transforms: dict[tuple[str, str], dt.Transform] = _load_reference_stucture_transforms(
                source_fs_model.reference_structure_transforms
            )
        else:
            reference_structure_transforms = _load_reference_stucture_transforms(
                fs_model.reference_structure_transforms
            )

        # Run the update
        _update(
            fs_model,
            datasets,
            reference_datasets,
            new_datasets,
            assemblies,
            xtalforms,
            dataset_assignments,
            ligand_neighbourhoods,
            alignability_graph,
            connected_components,
            ligand_neighbourhood_transforms,
            conformer_sites,
            conformer_site_transforms,
            canonical_sites,
            # canonical_site_transforms,
            xtalform_sites,
            reference_structure_transforms,
        )

    def process_all(self, option_json: str):
        options = Options.parse_file(option_json)

        # Initialize the output directory and create empty
        # jsons in it
        system_data = self.init(options.source_dir)

        # Add the datasources in the options json and add them to
        # the datasource json
        for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
            if datasource_type == "model_building":
                _add_model_building_dir_to_system_data(system_data, Path(datasource_dir))
            elif datasource_type == "manual":
                _add_manual_dir_to_system_data(system_data, Path(datasource_dir))

        # Add the PanDDAs in the options json and add them to the pandda json
        for pandda_dir in options.panddas:
            _add_pandda_to_system_data(system_data, Path(pandda_dir))

        # Copy the assembly json into the source directory (checking validity)
        # assemblies = Assemblies.read(Path(options.assemblies_json))
        # assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

        # Copy the xtalform json into the source directory (checking validity)
        xtalforms = XtalForms.read(Path(options.xtalforms_json))
        # xtalforms.save(Path(options.source_dir) / constants.XTALFORMS_FILE_NAME)

        # Parse the data sources and PanDDAs, matching ligands up to events
        system_data = _add_data_to_system_data(system_data)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

        # Build the alignment graph
        ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
            system_data,
            xtalforms,
            assigned_xtalforms,
        )

        num_neighbourhoods = len(ligand_neighbourhoods.ligand_neighbourhoods)
        logger.info(f"Found {num_neighbourhoods} ligand neighbourhoods")

        # Get alignability
        logger.info("Getting alignbaility matrix...!")
        alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)
        logger.info("Got alignability matrix!")

        # logger.debug(alignability_matrix)
        logger.debug("Alignability matrix shape: {alignability_matrix.shape}")

        # Generate the graph
        logger.info("Getting alignability graph...")
        g = get_graph(alignability_matrix, ligand_neighbourhoods)

        # Generate canonical, conformer and xtalform sites from the
        # alignment graph

        # Get the connected components
        logger.info("Getiting connected components...")
        connected_components = get_components(g)
        logger.info(f"Number of connected components: {len(connected_components)}")

        # Get the subsites from the connected components with overlap
        logger.info("Geting sites...")
        conformer_sites: ConformerSites = get_conformer_sites_from_components(
            connected_components, ligand_neighbourhoods
        )
        logger.info(f"Number of subsites: {len(conformer_sites.conformer_sites)}")

        # Merge the connected components with shared residues into sites
        logger.info("Getting sites...")
        _sites = get_sites_from_conformer_sites(conformer_sites, ligand_neighbourhoods)
        logger.info(f"Number of sites: {len(_sites)}")

        canonical_sites: CanonicalSites = CanonicalSites(
            site_ids=[s.id for s in _sites],
            sites=_sites,
            reference_site=_sites[0],
            reference_site_id=_sites[0].id,
        )

        # Get the xtalform sites
        # xtalform_sites = get_xtalform_sites_from_canonical_sites(
        #     canonical_sites,
        #     assigned_xtalforms,
        #     xtalforms,
        #     # assemblies,
        # )

        # Get the subsite transforms
        logger.info("Getting transfroms between subsites...")
        structures = get_structures(system_data)
        subsite_transforms = get_subsite_transforms(canonical_sites, structures)

        # Get the site transforms
        logger.info("Getting transforms between sites...")
        site_transforms = get_site_transforms(canonical_sites, structures)
        site_transforms = SiteTransforms(
            canonical_site_transform_ids=list(site_transforms.keys()),
            canonical_site_transforms=list(site_transforms.values()),
            conformer_site_transform_ids=list(subsite_transforms.keys()),
            conformer_site_transforms=list(subsite_transforms.values()),
        )
        # Fully specify the output now that the sites are known
        output = read_output(Path(options.source_dir))
        dataset_output_dict = {}
        for ligand_id in ligand_neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            if dtag not in dataset_output_dict:
                dataset_output = DatasetOutput(aligned_chain_output={})
                dataset_output_dict[dtag] = dataset_output
            else:
                dataset_output = dataset_output_dict[dtag]

            if chain not in dataset_output.aligned_chain_output:
                chain_output = ChainOutput(
                    aligned_ligands={},
                )
                dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
            else:
                chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

            chain_output.aligned_ligands[residue] = LigandOutput(
                aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

                chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

        # Save the output file
        output.dataset_output = dataset_output_dict

        # Align structures to each canonical site
        _align_structures_from_sites(
            structures,
            canonical_sites,
            conformer_sites,
            transforms,
            ligand_neighbourhoods,
            xtalforms,
            assigned_xtalforms,
            g,
            site_transforms,
            output,
        )

        # Align xmaps to each canonical site
        _align_xmaps(
            system_data,
            structures,
            canonical_sites,
            conformer_sites,
            ligand_neighbourhoods,
            g,
            transforms,
            site_transforms,
            output,
        )

    def process(self, option_json: str):
        options = Options.parse_file(option_json)

        # Initialize the output directory and create empty
        # jsons in it
        self.init(options.source_dir)

        # Add the datasources in the options json and add them to
        # the datasource json
        for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
            if datasource_type == "model_building":
                self.add_data_source(
                    options.source_dir,
                    datasource_dir,
                    source_type=datasource_type,
                )
            elif datasource_type == "manual":
                self.add_data_source(
                    options.source_dir,
                    datasource_dir,
                    source_type=datasource_type,
                )

        # Add the PanDDAs in the options json and add them to the pandda json
        for pandda_dir in options.panddas:
            self.add_pandda(options.source_dir, pandda_dir)

        # Copy the assembly json into the source directory (checking validity)
        assemblies = Assemblies.read(Path(options.assemblies_json))
        assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

        # Copy the xtalform json into the source directory (checking validity)
        xtalforms = XtalForms.read(Path(options.xtalforms_json))
        xtalforms.save(Path(options.source_dir) / constants.XTALFORMS_FILE_NAME)

        # Parse the data sources and PanDDAs, matching ligands up to events
        self.parse_data_sources(options.source_dir)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        self.assign_xtalforms(options.source_dir)

        # Build the alignment graph
        self.build_graph(options.source_dir)

        # Generate canonical, conformer and xtalform sites from the
        # alignment graph
        self.generate_sites_from_components(options.source_dir)

        # Fully specify the output now that the sites are known
        neighbourhoods = read_neighbourhoods(Path(options.source_dir))
        canonical_sites = CanonicalSites.read(Path(options.source_dir) / constants.CANONICAL_SITE_FILE)
        output = read_output(Path(options.source_dir))
        dataset_output_dict = {}
        for ligand_id in neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            if dtag not in dataset_output_dict:
                dataset_output = DatasetOutput(aligned_chain_output={})
                dataset_output_dict[dtag] = dataset_output
            else:
                dataset_output = dataset_output_dict[dtag]

            if chain not in dataset_output.aligned_chain_output:
                chain_output = ChainOutput(
                    aligned_ligands={},
                )
                dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
            else:
                chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

            chain_output.aligned_ligands[residue] = LigandOutput(
                aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

                chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                    output.aligned_dir
                    + "/"
                    + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                        dtag=dtag, chain=chain, residue=residue, site=site_id
                    )
                )

        # Save the output file
        output.dataset_output = dataset_output_dict
        save_output(output, Path(options.source_dir))

        # Align structures to each canonical site
        self.align_structures(options.source_dir)

        # Align xmaps to each canonical site
        self.align_xmaps(options.source_dir)

    def assign_xtalforms(self, source_dir: str):
        _source_dir = Path(source_dir)

        assemblies = Assemblies.read(_source_dir / constants.ASSEMBLIES_FILE_NAME)
        xtalforms = read_xtalforms(_source_dir)
        system_data = read_system_data(_source_dir)

        _assign_xtalforms(
            _source_dir,
            assemblies,
            xtalforms,
            system_data,  # dataset_xtalforms
        )

    def write_options_json(self, source_dir, options_json):
        _source_dir = Path(source_dir)
        _options_json = Path(options_json)

        system_data = read_system_data(_source_dir)
        options = Options(
            source_dir=source_dir,
            datasources=[ds.path for ds in system_data.datasources],
            panddas=[pandda.path for pandda in system_data.panddas],
        )
        with open(_options_json, "w") as f:
            f.write(options.json())

    def init(self, source_dir: str):
        _source_dir = Path(source_dir)

        if not _source_dir.exists():
            os.mkdir(_source_dir)

        system_data = SystemData(datasources=[], panddas=[], dataset_ids=[], datasets=[])

        save_data(system_data, _source_dir)

        output = Output(
            source_dir=str(_source_dir),
            system_data=str(constants.DATA_JSON_PATH),
            xtalforms=str(constants.XTALFORMS_FILE_NAME),
            assigned_xtalforms=str(constants.ASSIGNED_XTALFORMS_FILE_NAME),
            neighbourhoods=str(constants.NEIGHBOURHOODS_FILE_NAME),
            graph=str(constants.ALIGNABILITY_GRAPH_FILE_NAME),
            transforms=str(constants.TRANSFORMS_FILE_NAME),
            sites=str(constants.SITES_FILE_NAME),
            site_transforms=str(constants.SITES_TRANSFORMS_FILE_NAME),
            aligned_dir=str(constants.ALIGNED_STRUCTURES_DIR),
            dataset_output={},
        )
        if not (_source_dir / constants.ALIGNED_STRUCTURES_DIR).exists():
            os.mkdir(_source_dir / constants.ALIGNED_STRUCTURES_DIR)
        save_output(output, _source_dir)

        return system_data

    def add_data_source(
        self,
        source_dir: str,
        data_source_dir: str,
        source_type: str = "model_building",
    ):
        _source_dir = Path(source_dir)
        _data_source_dir = Path(data_source_dir)

        if source_type == "model_building":
            _add_model_building_dir(_source_dir, _data_source_dir)

        elif source_type == "manual":
            _add_manual_dir(_source_dir, _data_source_dir)

        else:
            raise Exception()

    def add_pandda(self, source_dir: str, pandda_dir: str):
        _source_dir = Path(source_dir)
        _pandda_dir = Path(pandda_dir)

        _add_pandda(_source_dir, _pandda_dir)

    def parse_data_sources(self, source_dir: str):
        _source_dir = Path(source_dir)

        _parse_data_sources(_source_dir)

    def open_site(
        self,
        option_json: str,
        site_id: int,
    ):
        options = Options.parse_file(option_json)
        output = Output.read(Path(options.source_dir) / constants.OUTPUT_JSON_PATH)

        _source_dir = Path(options.source_dir)
        script_path = _source_dir / "coot_script.py"
        script = ""
        script += 'if __name__ == "__main__": \n'
        script += '\tset_nomenclature_errors_on_read("ignore")\n'
        script += "\tset_recentre_on_read_pdb(0) \n"

        # str_dir = _source_dir / constants.ALIGNED_STRUCTURES_DIR

        # for site_dir in str_dir.glob("*"):
        #     if site_dir.name != f"{site_id}":
        #         continue

        #     for subsite_dir in site_dir.glob("*"):
        #         for pdb in subsite_dir.glob("*"):
        #             script += f'\tp = read_pdb("{pdb}")\n '
        #             script += cas_ligands()

        for dtag, dataset_output in output.dataset_output.items():
            for chain, chain_output in dataset_output.aligned_chain_output.items():
                for residue, residue_output in chain_output.aligned_ligands.items():
                    for _site_id, pdb in residue_output.aligned_structures.items():
                        logger.debug(_site_id)
                        if _site_id == site_id:
                            script += f'\tp = read_pdb("{options.source_dir}/{pdb}")\n '
                            script += cas_ligands()

        # for dataset in output.dataset_output

        logger.debug(script)

        with open(script_path, "w") as f:
            f.write(script)

        p = subprocess.Popen(f"coot --script {script_path}", shell=True)
        p.communicate()
        # os.remove(script_path)

    def merge_clusters(self, cluster_1: int, cluster_2: int, sites_path: str = "."):
        # sites = read_sites(sites_path)
        ...

    def suggest_merges(self):
        ...

    def pretty_print_dataset(self, source_dir: str):
        _source_dir = Path(source_dir)
        system_data = read_system_data(_source_dir)
        print(system_data)

    def align(
        self,
        system_data_dir: str,
        source_dir: str,
    ):
        self.build_system_data(system_data_dir, source_dir)
        self.build_graph(source_dir)
        self.generate_sites_from_components(source_dir)
        self.align_structures(source_dir)
        self.align_xmaps(source_dir)

    # def update(self, system_data_dir: str, source_dir: str):
    #
    #     # _source_dir: Path = Path(source_dir)
    #
    #     self.build_system_data(system_data_dir, source_dir)
    #     self.build_graph(source_dir)
    #     # self.update_sites(source_dir)
    #     self.align_structures(source_dir)
    #     self.align_xmaps(source_dir)

    def build_graph(
        self,
        source_dir: str,
    ):
        _source_dir: Path = Path(source_dir)

        build_alignment_graph(_source_dir)

    def build_system_data(self, system_data_dir: str, output_dir: str):
        print(output_dir)
        _system_data_dir: Path = Path(system_data_dir).resolve()
        _output_dir: Path = Path(output_dir).resolve()

        make_data_json_from_pandda_dir(_system_data_dir, _output_dir)

    def change_sites_reference(self, source_dir: str, site_id: int):
        _source_dir: Path = Path(source_dir)

        _change_sites_reference(_source_dir, site_id)

    def change_site_reference(self, source_dir: str, site_id: int, subsite_id: int):
        _source_dir: Path = Path(source_dir)

        _change_site_reference(_source_dir, site_id, subsite_id)

    def change_subsite_reference(
        self,
        source_dir: str,
        site_id: int,
        subsite_id: int,
        dtag: int,
        chain: str,
        residue: int,
    ):
        _source_dir: Path = Path(source_dir)

        _change_subsite_reference(_source_dir, site_id, subsite_id, dtag, chain, residue)

    def align_structures(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        assigned_xtalforms = read_assigned_xtalforms(_source_dir)
        canonical_sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
        conformer_sites: ConformerSites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)
        xtalforms = read_xtalforms(_source_dir)
        output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

        # get Structures
        structures = read_structures(system_data)

        # Align structures
        _align_structures_from_sites(
            structures,
            canonical_sites,
            conformer_sites,
            transforms,
            neighbourhoods,
            xtalforms,
            assigned_xtalforms,
            g,
            site_transforms,
            # _source_dir,
            output,
        )

    def align_xmaps(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        # xtalforms = XtalForms(xtalforms=[])
        sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
        conformer_sites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)
        output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

        # get Structures
        structures = read_structures(system_data)

        _align_xmaps(
            system_data,
            structures,
            sites,
            conformer_sites,
            neighbourhoods,
            g,
            transforms,
            site_transforms,
            output,
        )

    def generate_sites_from_components(self, source_dir: str):
        _source_dir: Path = Path(source_dir)

        _generate_sites_from_components(_source_dir)


if __name__ == "__main__":
    fire.Fire(CLI)
