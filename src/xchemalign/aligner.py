# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os
import traceback
import shutil
import time
from pathlib import Path

import yaml
import gemmi

from rich.traceback import install

# Local alignment imports
from ligand_neighbourhood_alignment import constants as lna_constants
from ligand_neighbourhood_alignment.align_xmaps import _align_xmaps

from ligand_neighbourhood_alignment import dt


from ligand_neighbourhood_alignment.cli import (
    _load_assemblies,
    _load_xtalforms,
    _load_dataset_assignments,
    _load_ligand_neighbourhoods,
    _load_alignability_graph,
    _load_connected_components,
    _load_ligand_neighbourhood_transforms,
    _load_conformer_sites,
    _load_conformer_site_transforms,
    _load_canonical_sites,
    _load_canonical_site_transforms,
    _load_xtalform_sites,
    _load_reference_stucture_transforms,
    _update,
)
from ligand_neighbourhood_alignment import alignment_heirarchy as ah

from xchemalign import utils
from xchemalign.utils import Constants
from xchemalign.pdb_xtal import PDBXtal

install(show_locals=True)


def try_make(path):
    if not Path(path).exists():
        os.mkdir(path)


def read_yaml(path):
    with open(path, 'r') as f:
        dic = yaml.safe_load(f)

    return dic


def path_to_relative_string(
    path,
    base_path,
):
    try:
        rel_path = path.relative_to(base_path)

        return str(rel_path)
    except AttributeError:
        return path
    except ValueError:
        return path


def traverse_dictionary(dic, func):
    for key, value in dic.items():
        try:
            traverse_dictionary(value, func)
        except AttributeError:
            dic[key] = func(value)
        except ValueError:
            dic[key] = func(value)


def _get_xmap_path_or_none(output_path, binding_event):
    xmap_file = binding_event.get(Constants.META_FILE)
    if xmap_file:
        return str(output_path / xmap_file)
    else:
        return None


def get_datasets_from_crystals(
    crystals,
    fs_model,
    output_path,
):
    # dataset_ids = [DatasetID(dtag=dtag) for dtag in crystals]
    # paths to files will be defined like this: upload_1/crystallographic_files/8dz1/8dz1.pdb
    # this is relative to the output_path variable that is defined from the metadata_file.yaml\
    datasets = {}
    reference_datasets = {}
    new_datasets = {}
    for dtag, crystal in crystals.items():
        mtz_file = crystal[Constants.META_XTAL_FILES].get(Constants.META_XTAL_MTZ, {}).get(Constants.META_FILE)
        if mtz_file is None:
            mtz_path = None
        else:
            mtz_path = str(output_path / mtz_file)

        dataset = dt.Dataset(
            dtag=dtag,
            pdb=str(output_path / crystal[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB][Constants.META_FILE]),
            xmap="",
            mtz=mtz_path,
            ligand_binding_events={
                (
                    str(dtag),
                    str(binding_event.get(Constants.META_PROT_CHAIN)),
                    str(binding_event.get(Constants.META_PROT_RES)),
                ): dt.LigandBindingEvent(
                    id=0,
                    dtag=str(dtag),
                    chain=str(binding_event.get(Constants.META_PROT_CHAIN)),
                    residue=str(binding_event.get(Constants.META_PROT_RES)),
                    xmap=_get_xmap_path_or_none(output_path, binding_event),
                )
                for binding_event in crystal[Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})
            },
        )
        datasets[dtag] = dataset
        if crystal[Constants.META_STATUS] == Constants.META_STATUS_NEW:
            new_datasets[dtag] = dataset
        if crystal[Constants.META_STATUS] == Constants.META_STATUS_SUPERSEDED:
            new_datasets[dtag] = dataset
        if dtag not in fs_model.alignments:
            new_datasets[dtag] = dataset
        if crystal.get(Constants.META_REFERENCE):
            reference_datasets[dtag] = dataset

    if (len(datasets) == 0) or (len(datasets) == 0):
        raise Exception

    for _dataset_id, _dataset in datasets.items():
        _num_binding_events = len(_dataset.ligand_binding_events)

    return datasets, reference_datasets, new_datasets


class Aligner:
    def __init__(self, dir, log_file=None, log_level=0):
        self.num_alignments = 0
        self.errors = []
        self.warnings = []

        self._find_version_dir(dir)  # sets self.working_dir and self.version_dir
        if not log_file:
            log_file = self.working_dir / 'upload-current' / 'aligner.log'
        self.logger = utils.Logger(logfile=log_file, level=log_level)
        self.logger.info("Using", self.version_dir, "as current version dir")
        self.base_dir = self.version_dir.parent  # e.g. path/to
        self.aligned_dir = (
            self.version_dir / Constants.META_ALIGNED_FILES
        )  # e.g. upload-current/upload_1/aligned_files
        self.xtal_dir = (
            self.version_dir / Constants.META_XTAL_FILES
        )  # e.g. upload-current/upload_1/crystallographic_files
        self.metadata_file = (
            self.version_dir / Constants.METADATA_COLLATOR_FILENAME
        )  # e.g. upload-current/upload_1/meta_collator.yaml
        self.assemblies_file = (
            self.base_dir / Constants.ASSEMBLIES_FILENAME
        )  # e.g. upload-current/upload_1/assemblies.yaml

    def _find_version_dir(self, dir):
        if dir:
            wd0 = Path(dir)
        else:
            wd0 = Path.cwd()

        wd1 = utils._verify_working_dir(wd0)
        if not wd1:
            self._log_error("Working dir " + wd0 + " is not valid")
            exit(1)

        self.working_dir = wd1
        current_dir = self.working_dir / 'upload-current'

        # check that we at least have an upload_1 dir
        if not (current_dir / 'upload_1').is_dir():
            self._log_error("Working dir " + dir + " does not contain and upload_? dirs")
            exit(1)

        # now find the latest upload_? dir
        i = 0
        while i < 100:
            i += 1
            version_dir = current_dir / ('upload_' + str(i))
            if version_dir.is_dir():
                self.version_dir = version_dir
            else:
                break

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def validate(self):
        if not self.version_dir.exists():
            self._log_error("version dir {} does not exist".format(self.version_dir))
        elif not self.version_dir.is_dir():
            self._log_error("version dir {} is not a directory".format(self.version_dir))
        else:
            output_meta_path = self.version_dir / Constants.METADATA_ALIGN_FILENAME
            if output_meta_path.exists():
                self._log_error(
                    "aligner output {} already exists. You must run aligner on clean output from collator.".format(
                        str(output_meta_path)
                    )
                )

            p = self.metadata_file
            if not p.exists():
                self._log_error("metadata file {} does not exist. Did the collator step run successfully?".format(p))
            if not p.is_file():
                self._log_error("metadata file {} is not a file. Did the collator step run successfully?".format(p))

            p = self.assemblies_file
            if not p.exists():
                self._log_error("assemblies.yaml file {} does not exist".format(p))
            elif not p.is_file():
                self._log_error("assemblies.yaml file {} is not a file".format(p))

        return len(self.errors), len(self.warnings)

    def run(self):
        self.logger.info("Running aligner...")
        input_meta = utils.read_config_file(str(self.metadata_file))
        collator_ver = input_meta.get(Constants.META_DATA_FORMAT_VERSION)
        if utils.DATA_FORMAT_VERSION != collator_ver:
            self._log_error(
                "Collator was run with a different data format version ({}). Current version is {}. "
                + "You must re-run collator so that it uses this version".format(
                    collator_ver, utils.DATA_FORMAT_VERSION
                )
            )
            exit(1)

        new_meta = self._perform_alignments(input_meta)
        self._write_output(input_meta, new_meta)

    def _write_output(self, collator_dict, aligner_dict):
        # keep a copy of the assemblies config
        self._copy_file_to_version_dir(self.assemblies_file)

        collator_dict[Constants.META_XTALFORMS] = aligner_dict[Constants.META_XTALFORMS]
        collator_dict[Constants.META_CONFORMER_SITES] = aligner_dict[Constants.META_CONFORMER_SITES]
        collator_dict[Constants.META_CANONICAL_SITES] = aligner_dict[Constants.META_CANONICAL_SITES]
        collator_dict[Constants.META_XTALFORM_SITES] = aligner_dict[Constants.META_XTALFORM_SITES]

        xtals = collator_dict[Constants.META_XTALS]
        # print(aligner_dict)
        for k, v in aligner_dict[Constants.META_XTALS].items():
            if k in xtals:
                xtals[k][Constants.META_ASSIGNED_XTALFORM] = v[Constants.META_ASSIGNED_XTALFORM]

            if Constants.META_ALIGNED_FILES in v:
                if k in xtals:
                    xtals[k][Constants.META_ALIGNED_FILES] = v[Constants.META_ALIGNED_FILES]
                    traverse_dictionary(
                        xtals[k][Constants.META_ALIGNED_FILES],
                        lambda x: path_to_relative_string(x, self.base_dir),
                    )

                else:
                    self.logger.warn('crystal {} not found in input. This is very strange.'.format(k))

        collator_dict[Constants.META_REFERENCE_ALIGNMENTS] = aligner_dict[Constants.META_REFERENCE_ALIGNMENTS]
        traverse_dictionary(
            collator_dict[Constants.META_REFERENCE_ALIGNMENTS],
            lambda x: path_to_relative_string(x, self.base_dir),
        )

        # remove this eventually
        with open(self.version_dir / 'aligner_tmp.yaml', "w") as stream:
            yaml.dump(aligner_dict, stream, sort_keys=False, default_flow_style=None)

        with open(self.version_dir / Constants.METADATA_ALIGN_FILENAME, "w") as stream:
            yaml.dump(collator_dict, stream, sort_keys=False, default_flow_style=None)

    def _copy_file_to_version_dir(self, file_path):
        f = shutil.copy2(file_path, self.version_dir)
        if not f:
            self.logger.warn("Failed to copy file {} to {}".format(file_path, self.version_dir))

    def _perform_alignments(self, meta):
        self.logger.info("Performing alignments")

        # Initialize the output directory and create empty
        # jsons in it

        # Add the datasources in the options json and add them to
        # the datasource json
        # visits = meta[lna_constants.META_INPUT]
        crystals = meta[Constants.META_XTALS]
        # Assert that
        if len(crystals) == 0:
            self.logger.error(f"Did not find any crystals in metadata file. Exiting.")
            raise Exception
        previous_version_dirs = meta.get(Constants.META_PREV_VERSION_DIRS)
        if len(previous_version_dirs) > 0:
            previous_output_path = self.base_dir / previous_version_dirs[-1]
        else:
            previous_output_path = None
        output_path = self.version_dir

        aligned_files_dir = output_path / Constants.META_ALIGNED_FILES
        if not aligned_files_dir.exists():
            os.mkdir(aligned_files_dir)

        # Load the previous output dir if there is one
        if previous_output_path:
            self.logger.info(f"Updating from previous directory: {previous_output_path}")
            source_fs_model = dt.FSModel.from_dir(previous_output_path)
        else:
            self.logger.info(f"First run! Not updating!")
            source_fs_model = None

        # Load the fs model for the new output dir
        fs_model = dt.FSModel.from_dir(output_path)
        fs_model.xtalforms = self.assemblies_file  # TODO change the name in LNA
        if source_fs_model:
            fs_model.alignments = source_fs_model.alignments
            fs_model.reference_alignments = source_fs_model.reference_alignments

        # Create output dir
        if not output_path.exists():
            os.mkdir(output_path)

        # Get the datasets
        datasets, reference_datasets, new_datasets = get_datasets_from_crystals(crystals, fs_model, self.base_dir)
        self.logger.info(f"Got {len(datasets)} datasets")
        for dtag, dataset in datasets.items():
            self.logger.info(f"Dataset {dtag} has {len(dataset.ligand_binding_events)} ligand binding events!")
        self.logger.info(f"Got {len(reference_datasets)} reference datasets")
        self.logger.info(f"Got {len(new_datasets)} new datasets")

        # Get assemblies
        if source_fs_model:
            assemblies: dict[str, dt.Assembly] = _load_assemblies(source_fs_model.xtalforms, self.assemblies_file)
        else:
            assemblies = _load_assemblies(fs_model.xtalforms, self.assemblies_file)
        self.logger.info(f"Got {len(assemblies)} assemblies")

        # # Get xtalforms
        if source_fs_model:
            xtalforms: dict[str, dt.XtalForm] = _load_xtalforms(source_fs_model.xtalforms, self.assemblies_file)
        else:
            xtalforms = _load_xtalforms(fs_model.xtalforms, self.assemblies_file)
        self.logger.info(f"Got {len(xtalforms)} xtalforms")

        # Get the dataset assignments
        if source_fs_model:
            dataset_assignments = _load_dataset_assignments(
                Path(source_fs_model.dataset_assignments), fail_if_not_found=True
            )
        else:
            dataset_assignments = _load_dataset_assignments(
                Path(fs_model.dataset_assignments), fail_if_not_found=False
            )

        # Get Ligand neighbourhoods
        if source_fs_model:
            ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = _load_ligand_neighbourhoods(
                source_fs_model.ligand_neighbourhoods, fail_if_not_found=True
            )
        else:
            ligand_neighbourhoods = _load_ligand_neighbourhoods(
                fs_model.ligand_neighbourhoods, fail_if_not_found=False
            )
        self.logger.info(f"Got {len(ligand_neighbourhoods)} ligand neighbourhoods")

        # Get alignability graph
        if source_fs_model:
            alignability_graph = _load_alignability_graph(source_fs_model.alignability_graph, fail_if_not_found=True)
        else:
            alignability_graph = _load_alignability_graph(fs_model.alignability_graph, fail_if_not_found=False)

        # Get the connected components
        if source_fs_model:
            connected_components = _load_connected_components(
                source_fs_model.connected_components, fail_if_not_found=True
            )
        else:
            connected_components = _load_connected_components(fs_model.connected_components, fail_if_not_found=False)

        if source_fs_model:
            self.logger.info(f"Have source fs model at {source_fs_model.ligand_neighbourhood_transforms}!")
            ligand_neighbourhood_transforms: dict[
                tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform
            ] = _load_ligand_neighbourhood_transforms(
                source_fs_model.ligand_neighbourhood_transforms, fail_if_not_found=True
            )
        else:
            ligand_neighbourhood_transforms = _load_ligand_neighbourhood_transforms(
                fs_model.ligand_neighbourhood_transforms, fail_if_not_found=False
            )
        # print(ligand_neighbourhood_transforms)

        # Get conformer sites
        if source_fs_model:
            conformer_sites: dict[str, dt.ConformerSite] = _load_conformer_sites(
                source_fs_model.conformer_sites, fail_if_not_found=True
            )
        else:
            conformer_sites = _load_conformer_sites(fs_model.conformer_sites, fail_if_not_found=False)

        #
        if source_fs_model:
            conformer_site_transforms: dict[tuple[str, str], dt.Transform] = _load_conformer_site_transforms(
                source_fs_model.conformer_site_transforms, fail_if_not_found=True
            )
        else:
            conformer_site_transforms = _load_conformer_site_transforms(
                fs_model.conformer_site_transforms, fail_if_not_found=False
            )

        # Get canonical sites
        if source_fs_model:
            canonical_sites: dict[str, dt.CanonicalSite] = _load_canonical_sites(
                source_fs_model.canonical_sites, fail_if_not_found=True
            )
        else:
            canonical_sites = _load_canonical_sites(fs_model.canonical_sites, fail_if_not_found=False)

        #
        if source_fs_model:
            canonical_site_transforms: dict[str, dt.Transform] = _load_canonical_site_transforms(
                source_fs_model.conformer_site_transforms, fail_if_not_found=True
            )
        else:
            canonical_site_transforms = _load_canonical_site_transforms(
                fs_model.conformer_site_transforms, fail_if_not_found=False
            )

        # Get xtalform sites
        if source_fs_model:
            xtalform_sites: dict[str, dt.XtalFormSite] = _load_xtalform_sites(
                source_fs_model.xtalform_sites, fail_if_not_found=True
            )
        else:
            xtalform_sites = _load_xtalform_sites(fs_model.xtalform_sites, fail_if_not_found=False)

        # Get reference structure transforms
        if source_fs_model:
            reference_structure_transforms: dict[tuple[str, str], dt.Transform] = _load_reference_stucture_transforms(
                source_fs_model.reference_structure_transforms, fail_if_not_found=True
            )
        else:
            reference_structure_transforms = _load_reference_stucture_transforms(
                fs_model.reference_structure_transforms, fail_if_not_found=False
            )

        # Get the assembly landmarks
        working_fs_model = fs_model
        if source_fs_model:
            working_fs_model = fs_model

        if working_fs_model.assembly_landmarks.exists():
            assembly_landmarks = ah.load_yaml(
                # CONOR_CHECK_THIS - self.base_dir was removed as a prefix to the first arg to be consistent
                # with other similar places, but no evidence that this is ever executed so not certain
                # that this is correct
                working_fs_model.assembly_landmarks,
                ah.dict_to_assembly_landmarks,
            )
        else:
            assembly_landmarks = {}

        # Get the assembly transforms
        if working_fs_model.assembly_transforms.exists():
            # CONOR_CHECK_THIS - self.base_dir was removed as a prefix to the first arg to be consistent
            # with other similar places, but no evidence that this is ever executed so not certain
            # that this is correct
            assembly_transforms = ah.load_yaml(working_fs_model.assembly_transforms, lambda x: x)
        else:
            assembly_transforms = {}

        # Run the update
        updated_fs_model = _update(
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
            xtalform_sites,
            reference_structure_transforms,
            assembly_landmarks,
            assembly_transforms,
            self.version_dir.name[7:],
        )

        # Update the metadata_file with aligned file locations and site information
        new_meta = {}

        # Add the xtalform information
        meta_xtalforms = {}
        xtalforms = read_yaml(updated_fs_model.xtalforms)
        for xtalform_id, xtalform in xtalforms[Constants.META_XTALFORMS].items():
            xtalform_reference = xtalform[Constants.META_REFERENCE]
            reference_structure = gemmi.read_structure(datasets[xtalform_reference].pdb)  # (xtalform_reference).pdb)
            reference_spacegroup = reference_structure.spacegroup_hm
            reference_unit_cell = reference_structure.cell

            meta_xtalforms[xtalform_id] = {
                Constants.META_XTALFORM_REFERENCE: xtalform_reference,
                Constants.META_XTALFORM_SPACEGROUP: reference_spacegroup,
                Constants.META_XTALFORM_CELL: {
                    "a": reference_unit_cell.a,
                    "b": reference_unit_cell.b,
                    "c": reference_unit_cell.c,
                    "alpha": reference_unit_cell.alpha,
                    "beta": reference_unit_cell.beta,
                    "gamma": reference_unit_cell.gamma,
                },
            }

        new_meta[Constants.META_XTALFORMS] = meta_xtalforms

        # Add the conformer sites
        conformer_sites = read_yaml(updated_fs_model.conformer_sites)
        new_meta[Constants.META_CONFORMER_SITES] = conformer_sites

        # Add the canonical sites
        canonical_sites = read_yaml(updated_fs_model.canonical_sites)
        new_meta[Constants.META_CANONICAL_SITES] = canonical_sites

        # Add the xtalform sites - note the chain is that of the original crystal structure, NOT the assembly
        xtalform_sites = read_yaml(updated_fs_model.xtalform_sites)
        new_meta[Constants.META_XTALFORM_SITES] = xtalform_sites

        # Add the output aligned files
        assigned_xtalforms = read_yaml(updated_fs_model.dataset_assignments)

        new_meta[Constants.META_XTALS] = {}
        for dtag, crystal in crystals.items():
            self.logger.info('looking at', dtag)

            new_meta[Constants.META_XTALS][dtag] = {}
            crystal_output = new_meta[Constants.META_XTALS][dtag]

            # Otherwise iterate the output data structure, adding the aligned structure,
            # artefacts, xmaps and event maps to the metadata_file
            assigned_xtalform = assigned_xtalforms[dtag]
            crystal_output[Constants.META_ASSIGNED_XTALFORM] = assigned_xtalform

            # Skip if no output for this dataset
            if dtag not in updated_fs_model.alignments:
                self.logger.warn('skipping {} as aligned structures not found'.format(dtag))
                continue

            # We capture the ligand binding events info as we need to know whether the event map file is present
            # This is a bit of a hack as the event map file location is generated by LNA even if there is no event map
            # so we need to know whether to actually include it in the metadata.
            # It would be better if LNA only included if it actually existed which would make the checking easier.
            event_map_dict_list = crystal.get(Constants.META_XTAL_FILES, {}).get(Constants.META_BINDING_EVENT, {})

            crystal_output[Constants.META_ALIGNED_FILES] = {}
            aligned_output = crystal_output[Constants.META_ALIGNED_FILES]
            dataset_output = updated_fs_model.alignments[dtag]
            # print(crystal)
            # print(dataset_output)
            # print(event_map_dict_list)
            for chain_name, chain_output in dataset_output.items():
                aligned_chain_output = aligned_output[chain_name] = {}
                i = 0
                for ligand_residue, ligand_output in chain_output.items():
                    aligned_ligand_output = aligned_chain_output[ligand_residue] = {}
                    for version, version_output in ligand_output.items():
                        aligned_version_output = aligned_ligand_output[version] = {}
                        for site_id, aligned_structure_path in version_output.aligned_structures.items():
                            # Is the event map file present?
                            # We do this assuming the order is the same as the info in the crystallographic files section
                            # But first do a check that event_map_dict_list contains the expected entry, as when the
                            # ligand residue number changes between versions then you get an error, so we try to
                            # handle this nicely. But there could be other things that cause this error.
                            if i >= len(event_map_dict_list):
                                msg = (
                                    'Unexpected number of event maps found for ligand '
                                    + ligand_residue
                                    + ' in crystal '
                                    + dtag
                                    + '. One cause for this is the ligand residue number changing between versions;'
                                    + ' if it\'s not that, the cause might be even more obscure.'
                                    + ' Talk to a developer, and good luck.'
                                    + ' If you can\'t correct this then you can add this crystal to the \`exclude\` list.'
                                )
                                self._log_error(msg)
                                exit(1)
                            event_map_present = True if Constants.META_FILE in event_map_dict_list[i] else False

                            aligned_artefacts_path = version_output.aligned_artefacts[site_id]
                            aligned_event_map_path = version_output.aligned_event_maps[site_id]
                            aligned_xmap_path = version_output.aligned_xmaps[site_id]
                            aligned_diff_map_path = version_output.aligned_diff_maps[site_id]

                            aligned_crystallographic_event_map_path = (
                                version_output.aligned_event_maps_crystallographic[site_id]
                            )
                            aligned_crystallographic_xmap_path = version_output.aligned_xmaps_crystallographic[site_id]
                            aligned_crystallographic_diff_map_path = version_output.aligned_diff_maps_crystallographic[
                                site_id
                            ]

                            aligned_version_output[site_id] = {
                                Constants.META_AIGNED_STRUCTURE: aligned_structure_path,
                                Constants.META_AIGNED_X_MAP: aligned_xmap_path,
                                Constants.META_AIGNED_DIFF_MAP: aligned_diff_map_path,
                                Constants.META_AIGNED_CRYSTALLOGRAPHIC_X_MAP: aligned_crystallographic_xmap_path,
                                Constants.META_AIGNED_CRYSTALLOGRAPHIC_DIFF_MAP: aligned_crystallographic_diff_map_path,
                            }
                            # if the event map is present then include it in the output
                            if event_map_present:
                                aligned_version_output[site_id][
                                    Constants.META_AIGNED_EVENT_MAP
                                ] = aligned_event_map_path
                                aligned_version_output[site_id][
                                    Constants.META_AIGNED_CRYSTALLOGRAPHIC_EVENT_MAP
                                ] = aligned_crystallographic_event_map_path
                    i += 1

        ## Add the reference alignments
        new_meta[Constants.META_REFERENCE_ALIGNMENTS] = {}
        for dtag, crystal in crystals.items():
            # Skip if no output for this dataset
            if dtag not in updated_fs_model.reference_alignments:
                continue

            crystal_output = new_meta[Constants.META_REFERENCE_ALIGNMENTS][dtag] = {}

            # for dtag, canonical_site_aligned_files in fs_model.reference_alignments.items():
            for canonical_site_id, aligned_files in updated_fs_model.reference_alignments[dtag].items():
                crystal_output[canonical_site_id] = aligned_files

        new_meta[Constants.META_TRANSFORMS] = {}

        ## Get the observation to conformer site transforms
        ligand_neighbourhood_transforms = read_yaml(updated_fs_model.ligand_neighbourhood_transforms)
        new_meta[Constants.META_TRANSFORMS][
            Constants.META_TRANSFORMS_OBSERVATION_TO_CONFORMER_SITES
        ] = ligand_neighbourhood_transforms

        ## Get the conformer site to canonical site transforms
        conformer_site_transforms = read_yaml(updated_fs_model.conformer_site_transforms)
        new_meta[Constants.META_TRANSFORMS][
            Constants.META_TRANSFORMS_CONFORMER_SITES_TO_CANON
        ] = conformer_site_transforms
        num_extract_errors = self._extract_components(crystals, new_meta)
        if num_extract_errors == 1:
            self.logger.warn(
                "there was a problem extracting components for 1 aligned structure. See above for details"
            )
        elif num_extract_errors > 1:
            self.logger.warn(
                "there were problems extracting components for",
                num_extract_errors,
                "aligned structures. See above for details",
            )

        # cleanup empty aligned files dirs
        empty_dir_count = 0
        for d in aligned_files_dir.iterdir():
            if len(list(d.iterdir())) == 0:
                empty_dir_count += 1
                d.rmdir()
        self.logger.info('removing {} empty aligned_files dirs'.format(empty_dir_count))
        return new_meta

    def _extract_components(self, crystals, aligner_meta):
        """
        Extract out the required forms of the molecules.
        1. *_apo.pdb - the aligned structure without the ligand
        2. *_apo_solv.pdb - the aligned solvent molecules only
        3. *_apo_desolv.pdb - the aligned structure protein chain only
        4. *_ligand.mol - molfile of the ligand
        5. *_ligand.pdb - PDB of the ligand

        :param meta:
        :return:
        """

        self.logger.info('extracting components')

        num_errors = 0
        num_pdbs = 0
        for k1, v1 in aligner_meta.get(Constants.META_XTALS, {}).items():  # k = xtal
            if Constants.META_ALIGNED_FILES in v1:
                self.logger.info('handling', k1)
                cif_file = (
                    crystals.get(k1)
                    .get(Constants.META_XTAL_FILES, {})
                    .get(Constants.META_XTAL_CIF, {})
                    .get(Constants.META_FILE)
                )
                for k2, v2 in v1[Constants.META_ALIGNED_FILES].items():  # chain
                    for k3, v3 in v2.items():  # ligand
                        for k4, v4 in v3.items():  # version
                            for k5, v5 in v4.items():  # conf site
                                pdb = v5[Constants.META_AIGNED_STRUCTURE]
                                self.logger.info("extracting components", k1, k2, k3, k4, k5, pdb)
                                # pth = self.version_dir / pdb
                                pth = Path(pdb)
                                if not pth.is_file():
                                    self.logger.error("can't find file", pth)
                                    num_errors += 1
                                else:
                                    pdbxtal = PDBXtal(pth, pth.parent, logger=self.logger)
                                    errs = pdbxtal.validate()
                                    if errs:
                                        self.logger.error("validation errors - can't extract components")
                                        num_errors += 1
                                    else:
                                        num_pdbs += 1
                                        pdbxtal.create_apo_file()
                                        pdbxtal.create_apo_solv_desolv()

                                        v5[Constants.META_PDB_APO] = str(pdbxtal.apo_file.relative_to(self.base_dir))
                                        v5[Constants.META_PDB_APO_SOLV] = str(
                                            pdbxtal.apo_solv_file.relative_to(self.base_dir)
                                        )
                                        v5[Constants.META_PDB_APO_DESOLV] = str(
                                            pdbxtal.apo_desolv_file.relative_to(self.base_dir)
                                        )
                                        if cif_file:
                                            try:
                                                pdbxtal.create_ligands(k2, k3, str(self.base_dir / cif_file))
                                                v5[Constants.META_LIGAND_MOL] = (
                                                    str(pdbxtal.ligand_base_file.relative_to(self.base_dir)) + '.mol'
                                                )
                                                v5[Constants.META_LIGAND_SDF] = (
                                                    str(pdbxtal.ligand_base_file.relative_to(self.base_dir)) + '.sdf'
                                                )
                                                v5[Constants.META_LIGAND_PDB] = (
                                                    str(pdbxtal.ligand_base_file.relative_to(self.base_dir)) + '.pdb'
                                                )
                                                v5[Constants.META_LIGAND_SMILES] = (
                                                    str(pdbxtal.ligand_base_file.relative_to(self.base_dir)) + '.smi'
                                                )
                                                v5[Constants.META_LIGAND_NAME] = pdbxtal.ligand_name
                                                v5[Constants.META_LIGAND_SMILES_STRING] = pdbxtal.smiles
                                            except:
                                                num_errors += 1
                                                self.logger.warn(
                                                    "failed to create ligand for",
                                                    k1,
                                                    "Check that the ligand in PDB file and the CIF file are compatible",
                                                )
                                                traceback.print_exc()

        self.num_alignments = num_pdbs
        return num_errors


def main():
    parser = argparse.ArgumentParser(description="aligner")

    parser.add_argument("-d", "--dir", help="Working directory")

    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()

    if args.dir:
        log = str(Path(args.dir) / 'aligner.log')
    else:
        log = 'aligner.log'

    try:
        a = Aligner(args.dir, log_file=log, log_level=args.log_level)
        logger = a.logger
        utils.LOG = logger

        num_errors, num_warnings = a.validate()

        if not args.validate:
            if num_errors:
                logger.error("There are errors, cannot continue")
                exit(1)
            else:
                t0 = time.time()
                a.run()
                t1 = time.time()
                logger.info("Handled {} alignments in {} secs".format(a.num_alignments, round(t1 - t0)))
                # write a summary of errors and warnings
                logger.report()
                logger.close()
                if logger.logfilename:
                    to_path = a.version_dir / 'aligner.log'
                    print("copying log file", logger.logfilename, "to", to_path)
                    f = shutil.copy2(logger.logfilename, to_path)
                    if not f:
                        print("Failed to copy log file {} to {}".format(logger.logfilename, to_path))
    except:
        # uncaught exception
        tb = traceback.format_exc()
        if args.dir:
            wd = args.dir
        else:
            wd = Path.cwd()
        logger.error(
            "Unexpected fatal error occurred when running aligner\n"
            + tb
            + "\nPlease send this information to the XCA developers:\nLog file: "
            + str(logger.logfilename)
            + "\nWorking dir location: "
            + str(wd)
        )


if __name__ == "__main__":
    main()
