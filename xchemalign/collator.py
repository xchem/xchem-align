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
from pathlib import Path
import shutil
import datetime
import re
import yaml


import gemmi
import numpy as np
import pandas as pd


from rdkit import Chem

from . import utils, dbreader
from .utils import Constants


def generate_xtal_dir(input_path: Path, xtal_name: str):
    """
    Generate the directory with the crystal data
    :param input_path:
    :param xtal_name:
    :return: The Path to the base dir for this crystal
    """
    xtal_dir = input_path / Constants.DEFAULT_MODEL_BUILDING_DIR / xtal_name
    return xtal_dir


def expand_file_path(filepath: Path, default="."):
    if filepath.is_absolute():
        return filepath
    else:
        return Path(default) / filepath


class Input:
    def __init__(
        self,
        base_path: Path,
        input_dir_path: Path,
        type: str,
        soakdb_file_path,
        panddas_event_file_paths: list[Path],
        reference=False,
        logger=None,
    ):
        self.base_path = base_path
        self.input_dir_path = input_dir_path
        self.type = type
        self.soakdb_file_path = soakdb_file_path
        self.panddas_event_file_paths = panddas_event_file_paths
        self.errors = []
        self.warnings = []
        self.reference = reference
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def get_input_dir_path(self, expand=True):
        if expand:
            return self.base_path / self.input_dir_path
        else:
            return self.input_dir_path

    def get_soakdb_file_path(self, expand=True):
        if expand:
            return self.base_path / self.input_dir_path / self.soakdb_file_path
        else:
            return self.soakdb_file_path

    def validate(self):
        if not self.base_path:
            self.errors.append("base path must be defined")
        elif not self.base_path.exists():
            self.errors.append("base_path does not exist: {}".format(self.base_path))
        elif not self.base_path.is_dir():
            self.errors.append("base_path is not a directory: {}".format(self.base_path))

        if not self.input_dir_path:
            self.errors.append("input_dir_path is not defined")
        elif self.input_dir_path.is_absolute():
            self.errors.append("input_path must be a path relative to base_path")
        else:
            p = self.get_input_dir_path()
            self.logger.info("testing", p)
            if not p.exists():
                self.errors.append("input_dir_path does not exist: {}".format(p))
            elif not p.is_dir():
                self.errors.append("input_dir_path is not a directory: {}".format(p))

        if self.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            if not self.soakdb_file_path:
                self.errors.append("soakdb_file_path is not defined")
            elif self.soakdb_file_path.is_absolute():
                self.errors.append("soakdb_file_path must be a path relative to input_path")
            else:
                p = self.get_soakdb_file_path()
                if not p.exists():
                    self.errors.append("soakdb_file_path does not exist: {}".format(p))
                elif not p.is_file():
                    self.errors.append("soakdb_file_path is not a file: {}".format(p))

        return len(self.errors), len(self.warnings)


class Collator:
    def __init__(self, config_file, logger=None):
        self.errors = []
        self.warnings = []
        self.config_file = config_file

        config = utils.read_config_file(config_file)
        self.config = config

        self.base_path = utils.find_path(config, Constants.CONFIG_BASE_DIR)
        self.output_path = utils.find_path(config, Constants.CONFIG_OUTPUT_DIR)

        self.target_name = utils.find_property(config, Constants.CONFIG_TARGET_NAME)

        self.version_number = None
        self.version_dir = None
        self.previous_version_dirs = []
        self.meta_history = []
        self.all_xtals = None
        self.new_or_updated_xtals = None
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

        self.inputs = []
        inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
        self.logger.info("found {} inputs".format(len(inputs)))
        if inputs:
            for input in inputs:
                input_path = utils.find_path(input, Constants.CONFIG_DIR)
                type = utils.find_property(input, Constants.CONFIG_TYPE)
                if type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                    soakdb_path = utils.find_path(
                        input, Constants.CONFIG_SOAKDB, default=Constants.DEFAULT_SOAKDB_PATH
                    )
                    panddas_csvs = utils.find_property(input, Constants.META_BINDING_EVENT)
                    if panddas_csvs:
                        panddas_paths = [Path(p) for p in panddas_csvs]
                    else:
                        panddas_paths = []

                    self.logger.info("adding input", input_path)
                    self.inputs.append(
                        Input(self.base_path, input_path, type, soakdb_path, panddas_paths, logger=self.logger)
                    )

                elif type == Constants.CONFIG_TYPE_MANUAL:
                    self.logger.info("adding input", input_path)
                    self.inputs.append(Input(self.base_path, input_path, type, None, [], logger=self.logger))
                else:
                    raise ValueError("unexpected input type:", type)

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def _log_warning(self, msg):
        self.logger.warn(msg)
        self.warnings.append(msg)

    def validate(self):
        v_dir = self.read_versions()
        if not v_dir:
            self.logger.error("Error with version dir. Please fix and try again.")
            return None
        self.logger.info("Using version dir {}".format(v_dir))

        self.logger.info("validating paths")
        num_errors, num_warnings = self.validate_paths()

        if num_errors == 0:
            self.logger.info("validating data")
            meta = self.validate_data()
        else:
            meta = {}

        return meta, len(self.errors), len(self.warnings)

    def validate_paths(self):
        self.errors.clear()
        self.warnings.clear()

        if not self.target_name:
            self._log_error("target_name not defined")
        else:
            if len(self.target_name) < 4:
                self._log_error("target_name must have at least 4 characters: " + self.target_name)
            else:
                x = re.search("^[A-Za-z]+[A-Za-z0-9_\\-]*$", self.target_name)
                if not x:
                    self._log_error("Invalid target_name: " + self.target_name)

        if self.output_path:
            if not self.output_path.exists():
                self._log_error("output_dir does not exist: {}".format(self.output_path))
            elif not self.output_path.is_dir():
                self._log_error("output_dir argument is not a directory: {}".format(self.output_path))
        else:
            self._log_error("Output dir not defined in the config file")

        num_input_errors = 0
        num_input_warnings = 0
        if self.inputs:
            for input in self.inputs:
                errors, warnings = input.validate()
                for error in input.errors:
                    self.logger.error(error)
                    num_input_errors += 1
                for warning in input.warnings:
                    self.logger.warn(warning)
                    num_input_warnings += 1

        return len(self.errors) + num_input_errors, len(self.warnings) + num_input_warnings

    def validate_data(self):
        """
        Read info from the SoakDB database and verify that the necessary entries are present
        :return: The generated metadata
        """
        crystals = {}
        input_dirs = []
        prev_version_dirs_str = [str(d) for d in self.previous_version_dirs]
        meta = {
            Constants.META_RUN_ON: str(datetime.datetime.now()),
            Constants.META_INPUT_DIRS: input_dirs,
            Constants.CONFIG_OUTPUT_DIR: str(self.output_path),
            Constants.META_VERSION_NUM: self.version_number,
            Constants.META_VERSION_DIR: str(self.version_dir),
            Constants.META_PREV_VERSION_DIRS: prev_version_dirs_str,
            Constants.META_XTALS: crystals,
        }

        for input in self.inputs:
            input_dirs.append(str(input.get_input_dir_path()))
            self._validate_input(input, crystals)

        return meta

    def _validate_input(self, input, crystals):
        if input.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            self._validate_soakdb_input(input, crystals)
        elif input.type == Constants.CONFIG_TYPE_MANUAL:
            self._validate_manual_input(input, crystals)
        else:
            raise ValueError("unexpected input type:", input.type)

    def _validate_soakdb_input(self, input, crystals):
        dbfile = input.get_soakdb_file_path()
        self.logger.info("opening DB file:", dbfile)
        df = dbreader.filter_dbmeta(dbfile)
        count = 0
        processed = 0
        num_pdb_files = 0
        num_mtz_files = 0
        num_cif_files = 0
        ref_datasets = set(self.config.get(Constants.CONFIG_REF_DATASETS, []))

        for index, row in df.iterrows():
            count += 1
            xtal_name = row[Constants.SOAKDB_XTAL_NAME]
            xtal_dir = generate_xtal_dir(input.input_dir_path, xtal_name)
            # print('xtal_dir:', xtal_dir)
            if not xtal_name:
                self._log_error("Crystal name not defined, cannot process row {}".format(xtal_name))
            else:
                # self.logger.info('Processing crystal {} {}'.format(count, xtal_name))
                missing_files = 0
                expanded_files = []
                colname = Constants.SOAKDB_COL_PDB
                file = row[colname]
                # RefinementPDB_latest file names are specified as absolute file names, but need to be handled as
                # relative to the base_path
                if not file:
                    expanded_files.append(None)
                    self._log_warning("PDB entry {} for {} not defined in SoakDB".format(colname, xtal_name))
                else:
                    # print('handling', colname, file)
                    inputpath = utils.make_path_relative(Path(file))
                    full_inputpath = self.base_path / inputpath
                    # print('generated', full_inputpath)
                    ok = full_inputpath.exists()
                    if ok:
                        num_pdb_files += 1
                        expanded_files.append(inputpath)
                    else:
                        expanded_files.append(None)
                        missing_files += 1
                        self._log_warning("PDB file for {} not found: {}".format(xtal_name, full_inputpath))

                    # if we have a PDB file then continue to look for the others
                    colname = Constants.SOAKDB_COL_MTZ
                    file = row[colname]
                    # RefinementMTZ_latest file names are specified as absolute file names, but need to be handled as
                    # relative to the base_path
                    if file:
                        inputpath = utils.make_path_relative(Path(file))
                        full_inputpath = self.base_path / inputpath
                        # print('generated', full_inputpath)
                        ok = full_inputpath.exists()
                        if ok:
                            num_mtz_files += 1
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            self._log_warning("MTZ file for {} not found: {}".format(xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning("MTZ entry {} for {} not defined in SoakDB".format(colname, xtal_name))

                    colname = Constants.SOAKDB_COL_CIF
                    file = row[colname]
                    # RefinementCIF file names are relative to the xtal_dir
                    if file:
                        p = Path(file)
                        if p.is_absolute():
                            inputpath = p.relative_to("/")
                        else:
                            inputpath = xtal_dir / file
                        full_inputpath = self.base_path / inputpath
                        # print('generated', full_inputpath)
                        ok = full_inputpath.exists()
                        if ok:
                            num_cif_files += 1
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            self._log_warning("CIF file for {} not found: {}".format(xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning("CIF entry {} for {} not defined in SoakDB".format(colname, xtal_name))

                    processed += 1
                    if xtal_name in crystals.keys():
                        self._log_warning("Crystal {} already exists, it's data will be overriden".format(xtal_name))

                    data = {}
                    if xtal_name in ref_datasets:
                        data[Constants.META_REFERENCE] = True
                    data[Constants.CONFIG_TYPE] = Constants.CONFIG_TYPE_MODEL_BUILDING
                    self.logger.info("adding crystal (model_building)", xtal_name)
                    crystals[xtal_name] = data
                    last_updated_date = row[Constants.SOAKDB_COL_LAST_UPDATED]
                    if last_updated_date:
                        dt_str = last_updated_date.strftime(utils._DATETIME_FORMAT)
                        data[Constants.META_LAST_UPDATED] = dt_str
                    digest = utils.gen_sha256(self.base_path / expanded_files[0])
                    f_data = {
                        Constants.META_XTAL_PDB: {
                            Constants.META_FILE: str(expanded_files[0]),
                            Constants.META_SHA256: digest,
                        }
                    }
                    if expanded_files[1]:
                        digest = utils.gen_sha256(self.base_path / expanded_files[1])
                        f_data[Constants.META_XTAL_MTZ] = {
                            Constants.META_FILE: str(expanded_files[1]),
                            Constants.META_SHA256: digest,
                        }
                    if expanded_files[2]:
                        digest = utils.gen_sha256(self.base_path / expanded_files[2])
                        f_data[Constants.META_XTAL_CIF] = {
                            Constants.META_FILE: str(expanded_files[2]),
                            Constants.META_SHA256: digest,
                        }
                    data[Constants.META_XTAL_FILES] = f_data

        print(data)

        self.logger.info("validator handled {} rows from database, {} were valid".format(count, processed))
        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding MTZ files".format(num_pdb_files, num_mtz_files)
            )
        if num_cif_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding CIF files".format(num_pdb_files, num_cif_files)
            )

    def _validate_manual_input(self, input, crystals):
        num_pdb_files = 0
        num_mtz_files = 0
        ref_datasets = set(self.config.get(Constants.CONFIG_REF_DATASETS, []))
        for child in (self.base_path / input.input_dir_path).iterdir():
            pdb = None
            mtz = None
            if child.is_dir():
                for file in child.iterdir():
                    if file.suffix == ".pdb":
                        pdb = file
                    if file.suffix == ".mtz":
                        mtz = file
                if pdb:
                    self.logger.info("adding crystal (manual)", child.name)
                    num_pdb_files += 1
                    digest = utils.gen_sha256(pdb)
                    data = {
                        Constants.META_XTAL_PDB: {
                            Constants.META_FILE: pdb.relative_to(self.base_path),
                            Constants.META_SHA256: digest,
                        }
                    }
                    if mtz:
                        digest = utils.gen_sha256(mtz)
                        data[Constants.META_XTAL_MTZ] = {
                            Constants.META_FILE: mtz.relative_to(self.base_path),
                            Constants.META_SHA256: digest,
                        }
                        num_mtz_files += 1
                    crystals[child.name] = {}
                    if child.name in ref_datasets:
                        crystals[child.name][Constants.META_REFERENCE] = True
                    crystals[child.name][Constants.CONFIG_TYPE] = Constants.CONFIG_TYPE_MANUAL
                    crystals[child.name][Constants.META_XTAL_FILES] = data

        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding MTZ files".format(num_pdb_files, num_mtz_files)
            )

    def read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = Constants.VERSION_DIR_PREFIX + str(version)
            v_path = self.output_path / v_dir
            if v_path.is_dir():
                version += 1
            else:
                break
        if version == 1:
            self.logger.error("No version directory found. Please create one named upload_1")
            return None

        # the working version dir is one less than the current value
        version -= 1
        self.logger.info("version is {}".format(version))
        v_dir = Constants.VERSION_DIR_PREFIX + str(version)
        v_path = self.output_path / v_dir

        # read the metadata from the earlier versions and record the version dirs
        if version > 1:
            for v in range(1, version):
                self.logger.info("reading metadata for version {}".format(v))
                dir_name = Constants.VERSION_DIR_PREFIX + str(v)
                meta = self.read_metadata(self.output_path / dir_name)
                self.meta_history.append(meta)
                self.previous_version_dirs.append(dir_name)

        self.logger.info("setting version dir to {}".format(v_dir))
        self.version_dir = Path(v_dir)
        self.version_number = version

        num_old_metas = len(self.meta_history)
        if num_old_metas:
            self.logger.info("found {} metadata files from previous versions".format(num_old_metas))

        return v_dir

    def read_metadata(self, version_dir):
        self.logger.info("reading metadata for version {}".format(version_dir))
        meta_file = version_dir / Constants.METADATA_XTAL_FILENAME

        meta = utils.read_config_file(str(meta_file))
        return meta

    def run(self, meta):
        self.logger.info("running collator...")

        self.logger.info("coping files ...")
        new_meta = self._copy_files(meta)
        self.logger.info("munging the history ...")
        all_xtals, new_xtals = self._munge_history(meta)
        self.logger.info("writing metadata ...")
        self._write_metadata(new_meta, all_xtals, new_xtals)
        self.logger.info("copying config ...")
        self._copy_config()
        self.logger.info("run complete")
        return new_meta

    def _copy_files(self, meta):
        cryst_path = self.version_dir / Constants.META_XTAL_FILES
        ext_cryst_path = self.output_path / cryst_path
        self.logger.info("using cryst_dir of", ext_cryst_path)
        if ext_cryst_path.exists():
            self.logger.info("removing old cryst_dir")
            shutil.rmtree(ext_cryst_path)
        self.logger.info("creating cryst_dir", ext_cryst_path)
        os.makedirs(ext_cryst_path)

        num_event_maps = 0

        event_tables = self._find_event_tables()

        for xtal_name, xtal in meta[Constants.META_XTALS].items():
            dir = cryst_path / xtal_name

            historical_xtal_data = self._collate_crystallographic_files_history(xtal_name)
            curr_xtal_data = xtal[Constants.META_XTAL_FILES]
            type = xtal[Constants.CONFIG_TYPE]
            files_to_copy = {}

            # handle the PDB file
            pdb = curr_xtal_data.get(Constants.META_XTAL_PDB)
            if pdb:
                pdb_input = self.base_path / pdb[Constants.META_FILE]
                if pdb_input.is_file():
                    digest = utils.gen_sha256(pdb_input)
                    old_digest = historical_xtal_data.get(Constants.META_XTAL_PDB, {}).get(Constants.META_SHA256)
                    if digest != old_digest:
                        # PDB is changed
                        pdb_name = xtal_name + ".pdb"
                        pdb_output = dir / pdb_name
                        files_to_copy[Constants.META_XTAL_PDB] = (pdb_input, pdb_output, digest)

                # handle the MTZ file
                mtz = curr_xtal_data.get(Constants.META_XTAL_MTZ)
                if mtz:
                    mtz_file = mtz[Constants.META_FILE]
                    mtz_input = self.base_path / mtz_file
                    if mtz_input.is_file():
                        digest = utils.gen_sha256(mtz_input)
                        old_digest = historical_xtal_data.get(Constants.META_XTAL_MTZ, {}).get(Constants.META_SHA256)
                        if digest != old_digest:
                            mtz_name = xtal_name + ".mtz"
                            mtz_output = dir / mtz_name
                            files_to_copy[Constants.META_XTAL_MTZ] = (mtz_input, mtz_output, digest)
                    elif type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                        self.logger.warn("mtz file {} not present".format(mtz_input))
                else:
                    self.logger.warn("MTZ entry missing for {}".format(xtal_name))

                # handle the CIF file
                cif = curr_xtal_data.get(Constants.META_XTAL_CIF)
                if cif:
                    cif_file = cif[Constants.META_FILE]
                    cif_input = self.base_path / cif_file
                    if cif_input.is_file():
                        digest = utils.gen_sha256(cif_input)
                        old_digest = historical_xtal_data.get(Constants.META_XTAL_CIF, {}).get(Constants.META_SHA256)
                        if digest != old_digest:
                            cif_name = xtal_name + ".cif"
                            cif_output = dir / cif_name
                            files_to_copy[Constants.META_XTAL_CIF] = (cif_input, cif_output, digest)
                elif type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                    self.logger.warn("CIF entry missing for {}".format(xtal_name))

                # handle panddas event maps
                hist_event_maps = {}
                for event_map_data in historical_xtal_data.get(Constants.META_BINDING_EVENT, []):
                    model = event_map_data.get(Constants.META_PROT_MODEL)
                    chain = event_map_data.get(Constants.META_PROT_CHAIN)
                    res = event_map_data.get(Constants.META_PROT_RES)
                    if model is not None and chain and res:
                        hist_event_maps[(model, chain, res)] = event_map_data

                best_event_map_paths = self.get_dataset_event_maps(xtal_name, pdb_input, event_tables)
                num_identical_historical_event_maps = 0
                event_maps_to_copy = {}
                if best_event_map_paths:
                    for k, tup in best_event_map_paths.items():
                        path = tup[0]
                        digest = utils.gen_sha256(path)
                        ccp4_output = cryst_path / xtal_name / "{}_{}_{}.ccp4".format(k[0], k[1], k[2])
                        event_maps_to_copy[k] = (path, ccp4_output, digest, k, tup[1], tup[2])
                        hist_data = hist_event_maps.get(k)
                        if hist_data:
                            if digest == hist_data.get(Constants.META_SHA256):
                                num_identical_historical_event_maps += 1

            else:
                self.logger.error("PDB entry missing for {}".format(xtal_name))
                return meta

            # now copy the files
            self.logger.info("{} has {} files to copy".format(xtal_name, len(files_to_copy)))
            fdata = files_to_copy.get(Constants.META_XTAL_PDB)
            data_to_add = {}
            if fdata:
                os.makedirs(self.output_path / dir)
                f = shutil.copy2(fdata[0], self.output_path / fdata[1], follow_symlinks=True)
                if not f:
                    self.logger.error("Failed to copy PDB file {} to {}".format(fdata[0], self.output_path / fdata[1]))
                else:
                    data_to_add[Constants.META_XTAL_PDB] = {
                        Constants.META_FILE: str(fdata[1]),
                        Constants.META_SHA256: fdata[2],
                    }
                    # copy MTZ file
                    fdata = files_to_copy.get(Constants.META_XTAL_MTZ)
                    if fdata:
                        f = shutil.copy2(fdata[0], self.output_path / fdata[1], follow_symlinks=True)
                        if not f:
                            self.logger.error(
                                "Failed to copy MTZ file {} to {}".format(fdata[0], self.output_path / fdata[1])
                            )
                        else:
                            data_to_add[Constants.META_XTAL_MTZ] = {
                                Constants.META_FILE: str(fdata[1]),
                                Constants.META_SHA256: fdata[2],
                            }
                    fdata = files_to_copy.get(Constants.META_XTAL_CIF)

                    # copy CIF file
                    if fdata:
                        f = shutil.copy2(fdata[0], self.output_path / fdata[1], follow_symlinks=True)
                        if not f:
                            self.logger.error(
                                "Failed to copy CIF file {} to {}".format(fdata[0], self.output_path / fdata[1])
                            )
                        else:
                            data_to_add[Constants.META_XTAL_CIF] = {
                                Constants.META_FILE: str(fdata[1]),
                                Constants.META_SHA256: fdata[2],
                            }
                            try:
                                mol = utils.gen_mol_from_cif(str(self.output_path / fdata[1]))
                                smi = Chem.MolToSmiles(mol)
                                data_to_add[Constants.META_XTAL_CIF][Constants.META_SMILES] = smi
                            except:
                                self.logger.warn('failed to generate SMILES for ligand {}'.format(xtal_name))

                            # copy event maps
                    if event_maps_to_copy:
                        if num_identical_historical_event_maps == len(event_maps_to_copy):
                            historical_data = historical_xtal_data.get(Constants.META_BINDING_EVENT)
                            data_to_add[Constants.META_BINDING_EVENT] = historical_data
                        else:
                            new_data = []
                            for key, to_copy in event_maps_to_copy.items():
                                # print('copying {} to {}'.format(to_copy[0], self.output_path / to_copy[1]))
                                f = shutil.copy2(to_copy[0], self.output_path / to_copy[1], follow_symlinks=True)
                                if not f:
                                    self.logger.error(
                                        "Failed to copy Panddas file {} to {}".format(
                                            to_copy[0], self.output_path / to_copy[1]
                                        )
                                    )
                                else:
                                    d = {
                                        Constants.META_FILE: str(to_copy[1]),
                                        Constants.META_SHA256: to_copy[2],
                                        Constants.META_PROT_MODEL: int(to_copy[3][0]),
                                        Constants.META_PROT_CHAIN: to_copy[3][1],
                                        Constants.META_PROT_RES: to_copy[3][2],
                                        Constants.META_PROT_INDEX: to_copy[4],
                                        Constants.META_PROT_BDC: to_copy[5],
                                    }
                                    new_data.append(d)

                            data_to_add[Constants.META_BINDING_EVENT] = new_data

            new_xtal_data = {}
            for k, v in historical_xtal_data.items():
                new_xtal_data[k] = v
            for k, v in data_to_add.items():
                new_xtal_data[k] = v
            xtal[Constants.META_XTAL_FILES] = new_xtal_data

        return meta

    def _find_event_tables(self):
        event_tables = {}
        for input in self.inputs:
            for panddas_file in input.panddas_event_file_paths:
                pfp = panddas_file
                # self.logger.info('Reading CSV:', pfp)
                df = pd.read_csv(input.get_input_dir_path() / pfp)
                event_tables[input.get_input_dir_path() / pfp] = df
        return event_tables

    def _find_pdb_in_history(self, xtal_name, xtal_data):
        for metad in reversed(self.meta_history):
            xtals = metad[Constants.META_XTALS]
            if xtal_name in xtals:
                data = xtals[xtal_name][Constants.META_XTAL_FILES]
                pdb = data[Constants.META_XTAL_PDB]
                sha256 = pdb[Constants.META_SHA256]
                # print('testing', xtal_data)
                if sha256 == xtal_data[Constants.META_XTAL_PDB][Constants.META_SHA256]:
                    return data
        return None

    def _collate_crystallographic_files_history(self, xtal_name):
        history = {}
        for metad in self.meta_history:
            xtals = metad[Constants.META_XTALS]
            if xtal_name in xtals:
                data = xtals[xtal_name][Constants.META_XTAL_FILES]
                for key, value in data.items():
                    history[key] = value

        return history

    def _munge_history(self, meta):
        all_xtals = {}
        new_or_updated_xtals = {}

        # get any user defined overrides
        overrides = self.config.get("overrides", {})

        count = 0
        for metad in self.meta_history:
            count += 1
            self.logger.info("munging metadata {}".format(count))
            xtals = metad[Constants.META_XTALS]
            total = 0
            for xtal_name, xtal_data in xtals.items():
                total += 1
                all_xtals[xtal_name] = xtal_data
            self.logger.info("metadata {} has {} items".format(count, total))

        count += 1
        self.logger.info("munging current metadata")
        xtals = meta["crystals"]
        total = 0
        for xtal_name, xtal_data in xtals.items():
            total += 1
            if xtal_name in all_xtals:
                old_xtal_data = all_xtals[xtal_name]
                status = self._get_xtal_status(xtal_name, old_xtal_data, xtal_data)
                xtal_data[Constants.META_STATUS] = status
                if status != Constants.META_STATUS_UNCHANGED:
                    new_or_updated_xtals[xtal_name] = xtal_data
            else:
                xtal_data[Constants.META_STATUS] = Constants.META_STATUS_NEW
                new_or_updated_xtals[xtal_name] = xtal_data
            all_xtals[xtal_name] = xtal_data

            # look for any deprecations
            xtals_overrides = overrides.get(Constants.META_XTALS, {})
            xtal_override = overrides.get(Constants.META_XTALS, {}).get(xtal_name)
            if xtal_override:
                status_override = xtal_override.get(Constants.META_STATUS)
                if status_override:
                    status = status_override.get(Constants.META_STATUS)
                    if status:
                        xtal_data[Constants.META_STATUS] = status
                        reason = status_override.get(Constants.META_REASON)
                        self.logger.info("status for xtal {} is overridden to be {}".format(xtal_name, status))
                        if reason:
                            xtal_data[Constants.META_REASON] = reason
                        else:
                            self.logger.warn("status is overridden, but no reason was given")
                    else:
                        self.logger.warn("status is declared to be overridden, but no new status was given")

        self.logger.info("metadata {} has {} items".format(count, total))
        self.logger.info(
            "munging resulted in {} total xtals, {} are new or updated".format(
                len(all_xtals), len(new_or_updated_xtals)
            )
        )

        self.all_xtals = all_xtals
        self.new_or_updated_xtals = new_or_updated_xtals
        return all_xtals, new_or_updated_xtals

    def _get_xtal_status(self, xtal_name, old_data, new_data):
        """
        Compare status using the SHA256 digest of the PDB file.
        This should be reliable, but will be subject to saying the entry is updated even if there is a trivial
        change.

        :param xtal_name:
        :param old_data:
        :param new_data:
        :return:
        """
        old_pdb_data = old_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        new_pdb_data = new_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        old_sha256 = old_pdb_data[Constants.META_SHA256]
        new_sha256 = new_pdb_data[Constants.META_SHA256]
        if old_sha256 == new_sha256:
            return Constants.META_STATUS_UNCHANGED
        else:
            return Constants.META_STATUS_SUPERSEDES

    def _get_xtal_status_old(self, xtal_name, old_data, new_data):
        """
        Compare status using the last_updated property.
        Even if this can be relied on, it is not present for manual entries.

        :param xtal_name:
        :param old_data:
        :param new_data:
        :return:
        """
        old_date = old_data.get(Constants.META_LAST_UPDATED)
        if old_date:
            new_date = new_data.get(Constants.META_LAST_UPDATED)
        if not old_date or not new_date:
            self.logger.warn("Dates not defined for {}, must assume xtal is updated".format(xtal_name))
            return Constants.META_STATUS_SUPERSEDES
        elif utils.to_datetime(new_date) > utils.to_datetime(old_date):
            self.logger.info("Xtal {} is updated".format(xtal_name))
            return Constants.META_STATUS_SUPERSEDES
        else:
            # self.logger.info('Xtal {} is unchanged'.format(xtal_name))
            return Constants.META_STATUS_UNCHANGED

    def _write_metadata(self, meta, all_xtals, new_xtals):
        f = self.output_path / self.version_dir / Constants.METADATA_XTAL_FILENAME
        with open(f, "w") as stream:
            yaml.dump(meta, stream, sort_keys=False)
        f = self.output_path / self.version_dir / "all_xtals.yaml"
        with open(f, "w") as stream:
            yaml.dump(all_xtals, stream, sort_keys=False)
            f = self.output_path / self.version_dir / "new_xtals.yaml"
        with open(f, "w") as stream:
            yaml.dump(new_xtals, stream, sort_keys=False)

    def _copy_config(self):
        f = shutil.copy2(self.config_file, self.output_path / self.version_dir)
        if not f:
            print("Failed to copy config file to {}".format((self.output_path / self.version_dir)))
            return False
        return True

    def get_ligand_coords(
        self,
        structure: gemmi.Structure,
    ) -> dict[tuple[str, str, str], np.array]:
        ligand_coords = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name in Constants.LIGAND_NAMES:
                        poss = []
                        for atom in residue:
                            pos = atom.pos
                            poss.append([pos.x, pos.y, pos.z])

                        arr = np.array(poss)
                        mean = np.mean(arr, axis=0)
                        ligand_coords[(model.name, chain.name, residue.seqid.num)] = mean

        return ligand_coords

    def get_closest_event_map(
        self,
        xtal_name: str,
        ligand_coord: np.array,
        event_tables: dict[Path, pd.DataFrame],
    ) -> Path:
        distances = {}
        data = {}
        for pandda_path, event_table in event_tables.items():
            # print('Processing', xtal_name, pandda_path)
            dataset_events = event_table[event_table[Constants.EVENT_TABLE_DTAG] == xtal_name]
            for idx, row in dataset_events.iterrows():
                event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
                bdc = row[Constants.EVENT_TABLE_BDC]
                x, y, z = row[Constants.EVENT_TABLE_X], row[Constants.EVENT_TABLE_Y], row[Constants.EVENT_TABLE_Z]
                distance = np.linalg.norm(np.array([x, y, z]).flatten() - ligand_coord.flatten())
                # print('Distance:', distance)
                for template in Constants.EVENT_MAP_TEMPLATES:
                    event_map_path = (
                        pandda_path.parent.parent
                        / Constants.PROCESSED_DATASETS_DIR
                        / xtal_name
                        / template.format(dtag=xtal_name, event_idx=event_idx, bdc=bdc)
                    )
                    if event_map_path.exists() and event_map_path.is_file():
                        distances[event_map_path] = distance
                        data[event_map_path] = (event_idx, bdc)
                        continue

        if distances:
            k = min(distances, key=lambda _key: distances[_key])
            return k, data[k]
        else:
            return None, None

    def get_dataset_event_maps(
        self, xtal_name: str, pdb_file: Path, event_tables: dict[Path, pd.DataFrame]
    ) -> dict[tuple[str, str, str], Path]:
        # Get the relevant structure
        # self.logger.info('Reading', xtal_name, pdb_file)
        # self.logger.info('Using {} event tables'.format(len(event_tables)))

        structure = gemmi.read_structure(str(pdb_file))

        # Get the coordinates of ligands
        ligand_coords = self.get_ligand_coords(structure)

        # Get the closest events within some reasonable radius
        closest_event_maps = {}
        for ligand_key, ligand_coord in ligand_coords.items():
            # print('coord:', ligand_coord)
            closest_event_map, data = self.get_closest_event_map(xtal_name, ligand_coord, event_tables)
            if closest_event_map:
                # print('closest:', closest_event_map)
                closest_event_maps[ligand_key] = (closest_event_map, data[0], data[1])

        return closest_event_maps


def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-c", "--config-file", default="config.yaml", help="Configuration file")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("-v", "--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("collator: ", args)

    c = Collator(args.config_file, logger=logger)

    meta, num_errors, num_warnings = c.validate()

    if not args.validate:
        if num_errors:
            print("There are errors, cannot continue")
            exit(1)
        else:
            c.run(meta)


if __name__ == "__main__":
    main()
