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
from src.ligand_neighbourhood_alignment.map_alignment import _align_xmaps

from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import _get_ligand_neighbourhood


# from ligand_neighbourhood_alignment.get_system_sites import get_system_sites
from src.ligand_neighbourhood_alignment.dep.build_alignment_graph import build_alignment_graph
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
from src.ligand_neighbourhood_alignment.structure_alignment import _align_structures_from_sites
from src.ligand_neighbourhood_alignment.dep.generate_sites_from_components import (  # get_xtalform_sites_from_canonical_sites,
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
from src.ligand_neighbourhood_alignment.dep.make_data_json import (
    get_ligand_binding_events_from_panddas,
    get_ligand_binding_events_from_structure,
    make_data_json_from_pandda_dir,
)
from src.ligand_neighbourhood_alignment.structure_alignment import _align_structure, _align_reference_structure
from src.ligand_neighbourhood_alignment.map_alignment import read_xmap, read_xmap_from_mtz, __align_xmap
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


# class CLI:
#     def schema(self, output_dir: str):
#         _output_dir = Path(output_dir)

#         if not _output_dir.exists():
#             os.mkdir(_output_dir)

#         save_schema(SystemData, _output_dir)
#         save_schema(LigandNeighbourhoods, _output_dir)
#         save_schema(CanonicalSites, _output_dir)
#         save_schema(Transforms, _output_dir)
#         save_schema(SiteTransforms, _output_dir)
#         save_schema(SystemData, _output_dir)
#         save_schema(SystemData, _output_dir)
#         save_schema(SystemData, _output_dir)

#     def update(self, options_json: str):
#         options = Options.parse_file(options_json)
#         logger.info(f"Input dir is: {options.source_dir}")
#         logger.info(f"Output dir is: {options.output_dir}")

#         if options.source_dir:
#             source_fs_model = dt.FSModel.from_dir(options.source_dir)
#         else:
#             source_fs_model = None

#         fs_model = dt.FSModel.from_dir(
#             options.output_dir,
#         )
#         if source_fs_model:
#             fs_model.alignments = source_fs_model.alignments
#             fs_model.reference_alignments = source_fs_model.reference_alignments

#         # Update the output fs model, creating flat symlinks to old data

#         if not Path(options.output_dir).exists():
#             os.mkdir(options.output_dir)

#         aligned_files_dir = Path(options.output_dir) / constants.ALIGNED_FILES_DIR
#         if not aligned_files_dir.exists():
#             os.mkdir(aligned_files_dir)

#         fs_model.symlink_old_data()

#         source_data_model = dt.SourceDataModel.from_fs_model(
#             fs_model, options.datasources, options.datasource_types, options.panddas
#         )

#         datasets, reference_datasets, new_datasets = source_data_model.get_datasets()

#         # Get assemblies
#         logger.info(f"Getting assemblies...")
#         if source_fs_model:
#             assemblies: dict[str, dt.Assembly] = _load_assemblies(
#                 source_fs_model.assemblies, Path(options.assemblies_json)
#             )
#         else:
#             assemblies = _load_assemblies(fs_model.assemblies, Path(options.assemblies_json))

#         # Get xtalforms
#         logger.info(f"Getting xtalforms...")
#         if source_fs_model:
#             xtalforms: dict[str, dt.XtalForm] = _load_xtalforms(
#                 source_fs_model.xtalforms, Path(options.xtalforms_json)
#             )
#         else:
#             xtalforms = _load_xtalforms(fs_model.xtalforms, Path(options.xtalforms_json))

#         # Get the dataset assignments
#         logger.info(f"Getting dataset assignments...")
#         if source_fs_model:
#             dataset_assignments = _load_dataset_assignments(Path(source_fs_model.dataset_assignments))
#         else:
#             dataset_assignments = _load_dataset_assignments(Path(fs_model.dataset_assignments))

#         # Get Ligand neighbourhoods
#         logger.info(f"Getting ligand neighbourhoods...")
#         if source_fs_model:
#             ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = _load_ligand_neighbourhoods(
#                 source_fs_model.ligand_neighbourhoods
#             )
#         else:
#             ligand_neighbourhoods = _load_ligand_neighbourhoods(fs_model.ligand_neighbourhoods)

#         # Get alignability graph
#         logger.info(f"Getting alignability graph...")
#         if source_fs_model:
#             alignability_graph = _load_alignability_graph(source_fs_model.alignability_graph)
#         else:
#             alignability_graph = _load_alignability_graph(fs_model.alignability_graph)

#         if source_fs_model:
#             connected_components = _load_connected_components(source_fs_model.connected_components)
#         else:
#             connected_components = _load_connected_components(fs_model.connected_components)

#         #
#         logger.info(f"Getting lighand neighbourhood transforms...")
#         if source_fs_model:
#             ligand_neighbourhood_transforms: dict[
#                 tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform
#             ] = _load_ligand_neighbourhood_transforms(source_fs_model.ligand_neighbourhood_transforms)
#         else:
#             ligand_neighbourhood_transforms = _load_ligand_neighbourhood_transforms(
#                 fs_model.ligand_neighbourhood_transforms
#             )

#         # Get conformer sites
#         logger.info(f"Getting conformer sites...")
#         if source_fs_model:
#             conformer_sites: dict[str, dt.ConformerSite] = _load_conformer_sites(source_fs_model.conformer_sites)
#         else:
#             conformer_sites = _load_conformer_sites(fs_model.conformer_sites)

#         #
#         logger.info(f"Getting conformer site transforms...")
#         if source_fs_model:
#             conformer_site_transforms: dict[tuple[str, str], dt.Transform] = _load_conformer_site_transforms(
#                 source_fs_model.conformer_site_transforms
#             )
#         else:
#             conformer_site_transforms = _load_conformer_site_transforms(fs_model.conformer_site_transforms)

#         # Get canonical sites
#         logger.info(f"Getting canonical sites...")
#         if source_fs_model:
#             canonical_sites: dict[str, dt.CanonicalSite] = _load_canonical_sites(source_fs_model.canonical_sites)
#         else:
#             canonical_sites = _load_canonical_sites(fs_model.canonical_sites)

#         # Get xtalform sites
#         logger.info(f"Getting xtalform sites...")
#         if source_fs_model:
#             xtalform_sites: dict[str, dt.XtalFormSite] = _load_xtalform_sites(source_fs_model.xtalform_sites)
#         else:
#             xtalform_sites = _load_xtalform_sites(fs_model.xtalform_sites)

#         # Get reference structure transforms
#         logger.info(f"Getting reference structure transforms...")
#         if source_fs_model:
#             reference_structure_transforms: dict[tuple[str, str], dt.Transform] = _load_reference_stucture_transforms(
#                 source_fs_model.reference_structure_transforms
#             )
#         else:
#             reference_structure_transforms = _load_reference_stucture_transforms(
#                 fs_model.reference_structure_transforms
#             )

#         # Run the update
#         _update(
#             fs_model,
#             datasets,
#             reference_datasets,
#             new_datasets,
#             assemblies,
#             xtalforms,
#             dataset_assignments,
#             ligand_neighbourhoods,
#             alignability_graph,
#             connected_components,
#             ligand_neighbourhood_transforms,
#             conformer_sites,
#             conformer_site_transforms,
#             canonical_sites,
#             # canonical_site_transforms,
#             xtalform_sites,
#             reference_structure_transforms,
#         )

#     def process_all(self, option_json: str):
#         options = Options.parse_file(option_json)

#         # Initialize the output directory and create empty
#         # jsons in it
#         system_data = self.init(options.source_dir)

#         # Add the datasources in the options json and add them to
#         # the datasource json
#         for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
#             if datasource_type == "model_building":
#                 _add_model_building_dir_to_system_data(system_data, Path(datasource_dir))
#             elif datasource_type == "manual":
#                 _add_manual_dir_to_system_data(system_data, Path(datasource_dir))

#         # Add the PanDDAs in the options json and add them to the pandda json
#         for pandda_dir in options.panddas:
#             _add_pandda_to_system_data(system_data, Path(pandda_dir))

#         # Copy the assembly json into the source directory (checking validity)

#         # Copy the xtalform json into the source directory (checking validity)
#         xtalforms = XtalForms.read(Path(options.xtalforms_json))

#         # Parse the data sources and PanDDAs, matching ligands up to events
#         system_data = _add_data_to_system_data(system_data)

#         # Assign each dataset to the clsoest xtalform and fail if this
#         # is not possible
#         assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

#         # Build the alignment graph
#         ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
#             system_data,
#             xtalforms,
#             assigned_xtalforms,
#         )

#         num_neighbourhoods = len(ligand_neighbourhoods.ligand_neighbourhoods)
#         logger.info(f"Found {num_neighbourhoods} ligand neighbourhoods")

#         # Get alignability
#         logger.info("Getting alignbaility matrix...!")
#         alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)
#         logger.info("Got alignability matrix!")

#         logger.debug("Alignability matrix shape: {alignability_matrix.shape}")

#         # Generate the graph
#         logger.info("Getting alignability graph...")
#         g = get_graph(alignability_matrix, ligand_neighbourhoods)

#         # Generate canonical, conformer and xtalform sites from the
#         # alignment graph

#         # Get the connected components
#         logger.info("Getiting connected components...")
#         connected_components = get_components(g)
#         logger.info(f"Number of connected components: {len(connected_components)}")

#         # Get the subsites from the connected components with overlap
#         logger.info("Geting sites...")
#         conformer_sites: ConformerSites = get_conformer_sites_from_components(
#             connected_components, ligand_neighbourhoods
#         )
#         logger.info(f"Number of subsites: {len(conformer_sites.conformer_sites)}")

#         # Merge the connected components with shared residues into sites
#         logger.info("Getting sites...")
#         _sites = get_sites_from_conformer_sites(conformer_sites, ligand_neighbourhoods)
#         logger.info(f"Number of sites: {len(_sites)}")

#         canonical_sites: CanonicalSites = CanonicalSites(
#             site_ids=[s.id for s in _sites],
#             sites=_sites,
#             reference_site=_sites[0],
#             reference_site_id=_sites[0].id,
#         )

#         # Get the subsite transforms
#         logger.info("Getting transfroms between subsites...")
#         structures = get_structures(system_data)
#         subsite_transforms = get_subsite_transforms(canonical_sites, structures)

#         # Get the site transforms
#         logger.info("Getting transforms between sites...")
#         site_transforms = get_site_transforms(canonical_sites, structures)
#         site_transforms = SiteTransforms(
#             canonical_site_transform_ids=list(site_transforms.keys()),
#             canonical_site_transforms=list(site_transforms.values()),
#             conformer_site_transform_ids=list(subsite_transforms.keys()),
#             conformer_site_transforms=list(subsite_transforms.values()),
#         )
#         # Fully specify the output now that the sites are known
#         output = read_output(Path(options.source_dir))
#         dataset_output_dict = {}
#         for ligand_id in ligand_neighbourhoods.ligand_ids:
#             dtag, chain, residue = (
#                 ligand_id.dtag,
#                 ligand_id.chain,
#                 ligand_id.residue,
#             )

#             if dtag not in dataset_output_dict:
#                 dataset_output = DatasetOutput(aligned_chain_output={})
#                 dataset_output_dict[dtag] = dataset_output
#             else:
#                 dataset_output = dataset_output_dict[dtag]

#             if chain not in dataset_output.aligned_chain_output:
#                 chain_output = ChainOutput(
#                     aligned_ligands={},
#                 )
#                 dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
#             else:
#                 chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

#             chain_output.aligned_ligands[residue] = LigandOutput(
#                 aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
#             )

#             # Add output for each canonical site that the ligand is aligned to
#             for site_id, site in canonical_sites.iter():
#                 if ligand_id not in site.members:
#                     continue

#                 chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#                 chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#                 chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
#                 )

#                 chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#         # Save the output file
#         output.dataset_output = dataset_output_dict

#         # Align structures to each canonical site
#         _align_structures_from_sites(
#             structures,
#             canonical_sites,
#             conformer_sites,
#             transforms,
#             ligand_neighbourhoods,
#             xtalforms,
#             assigned_xtalforms,
#             g,
#             site_transforms,
#             output,
#         )

#         # Align xmaps to each canonical site
#         _align_xmaps(
#             system_data,
#             structures,
#             canonical_sites,
#             conformer_sites,
#             ligand_neighbourhoods,
#             g,
#             transforms,
#             site_transforms,
#             output,
#         )

#     def process(self, option_json: str):
#         options = Options.parse_file(option_json)

#         # Initialize the output directory and create empty
#         # jsons in it
#         self.init(options.source_dir)

#         # Add the datasources in the options json and add them to
#         # the datasource json
#         for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
#             if datasource_type == "model_building":
#                 self.add_data_source(
#                     options.source_dir,
#                     datasource_dir,
#                     source_type=datasource_type,
#                 )
#             elif datasource_type == "manual":
#                 self.add_data_source(
#                     options.source_dir,
#                     datasource_dir,
#                     source_type=datasource_type,
#                 )

#         # Add the PanDDAs in the options json and add them to the pandda json
#         for pandda_dir in options.panddas:
#             self.add_pandda(options.source_dir, pandda_dir)

#         # Copy the assembly json into the source directory (checking validity)
#         assemblies = Assemblies.read(Path(options.assemblies_json))
#         assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

#         # Copy the xtalform json into the source directory (checking validity)
#         xtalforms = XtalForms.read(Path(options.xtalforms_json))
#         xtalforms.save(Path(options.source_dir) / constants.XTALFORMS_FILE_NAME)

#         # Parse the data sources and PanDDAs, matching ligands up to events
#         self.parse_data_sources(options.source_dir)

#         # Assign each dataset to the clsoest xtalform and fail if this
#         # is not possible
#         self.assign_xtalforms(options.source_dir)

#         # Build the alignment graph
#         self.build_graph(options.source_dir)

#         # Generate canonical, conformer and xtalform sites from the
#         # alignment graph
#         self.generate_sites_from_components(options.source_dir)

#         # Fully specify the output now that the sites are known
#         neighbourhoods = read_neighbourhoods(Path(options.source_dir))
#         canonical_sites = CanonicalSites.read(Path(options.source_dir) / constants.CANONICAL_SITE_FILE)
#         output = read_output(Path(options.source_dir))
#         dataset_output_dict = {}
#         for ligand_id in neighbourhoods.ligand_ids:
#             dtag, chain, residue = (
#                 ligand_id.dtag,
#                 ligand_id.chain,
#                 ligand_id.residue,
#             )

#             if dtag not in dataset_output_dict:
#                 dataset_output = DatasetOutput(aligned_chain_output={})
#                 dataset_output_dict[dtag] = dataset_output
#             else:
#                 dataset_output = dataset_output_dict[dtag]

#             if chain not in dataset_output.aligned_chain_output:
#                 chain_output = ChainOutput(
#                     aligned_ligands={},
#                 )
#                 dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
#             else:
#                 chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

#             chain_output.aligned_ligands[residue] = LigandOutput(
#                 aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
#             )

#             # Add output for each canonical site that the ligand is aligned to
#             for site_id, site in canonical_sites.iter():
#                 if ligand_id not in site.members:
#                     continue

#                 chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#                 chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#                 chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
#                 )

#                 chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
#                     output.aligned_dir
#                     + "/"
#                     + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
#                         dtag=dtag, chain=chain, residue=residue, site=site_id
#                     )
#                 )

#         # Save the output file
#         output.dataset_output = dataset_output_dict
#         save_output(output, Path(options.source_dir))

#         # Align structures to each canonical site
#         self.align_structures(options.source_dir)

#         # Align xmaps to each canonical site
#         self.align_xmaps(options.source_dir)

#     def assign_xtalforms(self, source_dir: str):
#         _source_dir = Path(source_dir)

#         assemblies = Assemblies.read(_source_dir / constants.ASSEMBLIES_FILE_NAME)
#         xtalforms = read_xtalforms(_source_dir)
#         system_data = read_system_data(_source_dir)

#         _assign_xtalforms(
#             _source_dir,
#             assemblies,
#             xtalforms,
#             system_data,  # dataset_xtalforms
#         )

#     def write_options_json(self, source_dir, options_json):
#         _source_dir = Path(source_dir)
#         _options_json = Path(options_json)

#         system_data = read_system_data(_source_dir)
#         options = Options(
#             source_dir=source_dir,
#             datasources=[ds.path for ds in system_data.datasources],
#             panddas=[pandda.path for pandda in system_data.panddas],
#         )
#         with open(_options_json, "w") as f:
#             f.write(options.json())

#     def init(self, source_dir: str):
#         _source_dir = Path(source_dir)

#         if not _source_dir.exists():
#             os.mkdir(_source_dir)

#         system_data = SystemData(datasources=[], panddas=[], dataset_ids=[], datasets=[])

#         save_data(system_data, _source_dir)

#         output = Output(
#             source_dir=str(_source_dir),
#             system_data=str(constants.DATA_JSON_PATH),
#             xtalforms=str(constants.XTALFORMS_FILE_NAME),
#             assigned_xtalforms=str(constants.ASSIGNED_XTALFORMS_FILE_NAME),
#             neighbourhoods=str(constants.NEIGHBOURHOODS_FILE_NAME),
#             graph=str(constants.ALIGNABILITY_GRAPH_FILE_NAME),
#             transforms=str(constants.TRANSFORMS_FILE_NAME),
#             sites=str(constants.SITES_FILE_NAME),
#             site_transforms=str(constants.SITES_TRANSFORMS_FILE_NAME),
#             aligned_dir=str(constants.ALIGNED_STRUCTURES_DIR),
#             dataset_output={},
#         )
#         if not (_source_dir / constants.ALIGNED_STRUCTURES_DIR).exists():
#             os.mkdir(_source_dir / constants.ALIGNED_STRUCTURES_DIR)
#         save_output(output, _source_dir)

#         return system_data

#     def add_data_source(
#         self,
#         source_dir: str,
#         data_source_dir: str,
#         source_type: str = "model_building",
#     ):
#         _source_dir = Path(source_dir)
#         _data_source_dir = Path(data_source_dir)

#         if source_type == "model_building":
#             _add_model_building_dir(_source_dir, _data_source_dir)

#         elif source_type == "manual":
#             _add_manual_dir(_source_dir, _data_source_dir)

#         else:
#             raise Exception()

#     def add_pandda(self, source_dir: str, pandda_dir: str):
#         _source_dir = Path(source_dir)
#         _pandda_dir = Path(pandda_dir)

#         _add_pandda(_source_dir, _pandda_dir)

#     def parse_data_sources(self, source_dir: str):
#         _source_dir = Path(source_dir)

#         _parse_data_sources(_source_dir)

#     def open_site(
#         self,
#         option_json: str,
#         site_id: int,
#     ):
#         options = Options.parse_file(option_json)
#         output = Output.read(Path(options.source_dir) / constants.OUTPUT_JSON_PATH)

#         _source_dir = Path(options.source_dir)
#         script_path = _source_dir / "coot_script.py"
#         script = ""
#         script += 'if __name__ == "__main__": \n'
#         script += '\tset_nomenclature_errors_on_read("ignore")\n'
#         script += "\tset_recentre_on_read_pdb(0) \n"

#         for dtag, dataset_output in output.dataset_output.items():
#             for chain, chain_output in dataset_output.aligned_chain_output.items():
#                 for residue, residue_output in chain_output.aligned_ligands.items():
#                     for _site_id, pdb in residue_output.aligned_structures.items():
#                         logger.debug(_site_id)
#                         if _site_id == site_id:
#                             script += f'\tp = read_pdb("{options.source_dir}/{pdb}")\n '
#                             script += cas_ligands()

#         # for dataset in output.dataset_output

#         logger.debug(script)

#         with open(script_path, "w") as f:
#             f.write(script)

#         p = subprocess.Popen(f"coot --script {script_path}", shell=True)
#         p.communicate()

#     def merge_clusters(self, cluster_1: int, cluster_2: int, sites_path: str = "."):
#         ...

#     def suggest_merges(self):
#         ...

#     def pretty_print_dataset(self, source_dir: str):
#         _source_dir = Path(source_dir)
#         system_data = read_system_data(_source_dir)
#         print(system_data)

#     def align(
#         self,
#         system_data_dir: str,
#         source_dir: str,
#     ):
#         self.build_system_data(system_data_dir, source_dir)
#         self.build_graph(source_dir)
#         self.generate_sites_from_components(source_dir)
#         self.align_structures(source_dir)
#         self.align_xmaps(source_dir)

#     def build_graph(
#         self,
#         source_dir: str,
#     ):
#         _source_dir: Path = Path(source_dir)

#         build_alignment_graph(_source_dir)

#     def build_system_data(self, system_data_dir: str, output_dir: str):
#         print(output_dir)
#         _system_data_dir: Path = Path(system_data_dir).resolve()
#         _output_dir: Path = Path(output_dir).resolve()

#         make_data_json_from_pandda_dir(_system_data_dir, _output_dir)

#     def change_sites_reference(self, source_dir: str, site_id: int):
#         _source_dir: Path = Path(source_dir)

#         _change_sites_reference(_source_dir, site_id)

#     def change_site_reference(self, source_dir: str, site_id: int, subsite_id: int):
#         _source_dir: Path = Path(source_dir)

#         _change_site_reference(_source_dir, site_id, subsite_id)

#     def change_subsite_reference(
#         self,
#         source_dir: str,
#         site_id: int,
#         subsite_id: int,
#         dtag: int,
#         chain: str,
#         residue: int,
#     ):
#         _source_dir: Path = Path(source_dir)

#         _change_subsite_reference(_source_dir, site_id, subsite_id, dtag, chain, residue)

#     def align_structures(self, source_dir: str):
#         _source_dir: Path = Path(source_dir)

#         g = read_graph(_source_dir)
#         transforms: Transforms = read_transforms(_source_dir)
#         neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
#         assigned_xtalforms = read_assigned_xtalforms(_source_dir)
#         canonical_sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
#         conformer_sites: ConformerSites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
#         system_data: SystemData = read_system_data(_source_dir)
#         site_transforms = read_site_transforms(_source_dir)
#         xtalforms = read_xtalforms(_source_dir)
#         output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

#         # get Structures
#         structures = read_structures(system_data)

#         # Align structures
#         _align_structures_from_sites(
#             structures,
#             canonical_sites,
#             conformer_sites,
#             transforms,
#             neighbourhoods,
#             xtalforms,
#             assigned_xtalforms,
#             g,
#             site_transforms,
#             # _source_dir,
#             output,
#         )

#     def align_xmaps(self, source_dir: str):
#         _source_dir: Path = Path(source_dir)

#         g = read_graph(_source_dir)
#         transforms: Transforms = read_transforms(_source_dir)
#         neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
#         sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
#         conformer_sites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
#         system_data: SystemData = read_system_data(_source_dir)
#         site_transforms = read_site_transforms(_source_dir)
#         output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

#         # get Structures
#         structures = read_structures(system_data)

#         _align_xmaps(
#             system_data,
#             structures,
#             sites,
#             conformer_sites,
#             neighbourhoods,
#             g,
#             transforms,
#             site_transforms,
#             output,
#         )

#     def generate_sites_from_components(self, source_dir: str):
#         _source_dir: Path = Path(source_dir)

#         _generate_sites_from_components(_source_dir)


# if __name__ == "__main__":
#     fire.Fire(CLI)
