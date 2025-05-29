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
import math
import re
import time
import traceback
from distutils import dir_util
from distutils.dir_util import copy_tree

import yaml

import gemmi
import numpy as np
import pandas as pd

from rdkit import Chem

from xchemalign import dbreader
from xchemalign import utils
from xchemalign import repo_info
from xchemalign import setup
from xchemalign.utils import Constants


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
        exclude: list[str],
        code_prefix=None,
        code_prefix_tooltip=None,
        reference=False,
        logger=None,
    ):
        self.base_path = base_path
        self.input_dir_path = input_dir_path
        self.type = type
        self.soakdb_file_path = soakdb_file_path
        self.panddas_event_file_paths = panddas_event_file_paths
        self.exclude = exclude
        self.errors = []
        self.warnings = []
        self.code_prefix = code_prefix
        self.code_prefix_tooltip = code_prefix_tooltip
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
            self.errors.append(
                "inputs.dir must be a path relative to base_path. What was specified was", self.input_dir_path
            )
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

        if self.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            if self.code_prefix is None:
                self.errors.append("code_prefix property is not defined")
            if self.code_prefix_tooltip is None:
                self.errors.append("code_prefix_tooltip property is not defined")

        return len(self.errors), len(self.warnings)


class Collator:
    def __init__(self, working_dir, log_file=None, log_level=0, include_git_info=False):
        self.errors = []
        self.warnings = []

        self.working_dir = Path(working_dir)
        self.config_file = self.working_dir / "upload-current" / "config.yaml"
        if not self.config_file.is_file():
            print(self.config_file, "not found")
        self.include_git_info = include_git_info

        config = utils.read_config_file(self.config_file)
        self.config = config

        self.base_path = utils.find_path(config, Constants.CONFIG_BASE_DIR)
        self.output_path = self.working_dir / "upload-current"

        self.target_name = utils.find_property(config, Constants.CONFIG_TARGET_NAME)

        self.panddas_missing_ok = utils.find_property(config, Constants.META_PANDDAS_MISSING_OK, default=[])

        self.version_number = None
        self.version_dir = None
        self.previous_version_dirs = []
        self.meta_history = []
        self.all_xtals = None
        self.rejected_xtals = set()
        self.all_excludes = set()
        self.new_or_updated_xtals = None
        self.compound_codes = {}
        self.compound_smiles = {}
        self.num_crystals = 0

        if not log_file:
            log_file = self.output_path / "collator.log"
        self.logger = utils.Logger(logfile=log_file, level=log_level)
        self.log_file = log_file

        self.inputs = []
        inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
        self.logger.info("found {} inputs".format(len(inputs)))

        if inputs:
            for input in inputs:
                # Determine which datasets to exclude
                excluded_datasets = utils.find_property(input, Constants.CONFIG_EXCLUDE)
                if excluded_datasets:
                    for excl in excluded_datasets:
                        self.all_excludes.add(excl)
                else:
                    excluded_datasets = []

                self.logger.info("Excluded datasets:", excluded_datasets)

                input_path = utils.find_path(input, Constants.CONFIG_DIR)
                type = utils.find_property(input, Constants.CONFIG_TYPE)
                code_prefix = utils.find_property(input, Constants.CONFIG_CODE_PREFIX)
                code_prefix_tooltip = utils.find_property(input, Constants.CONFIG_CODE_PREFIX_TOOLTIP)
                if type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                    soakdb_path = utils.find_path(
                        input, Constants.CONFIG_SOAKDB, default=Constants.DEFAULT_SOAKDB_PATH
                    )
                    panddas_csvs = utils.find_property(input, Constants.CONFIG_PANDDAS_EVENT_FILES)
                    if panddas_csvs:
                        panddas_paths = [Path(p) for p in panddas_csvs]
                    else:
                        panddas_paths = []

                    self.logger.info("adding input", input_path, "with", len(panddas_paths), "panddas event maps")
                    self.inputs.append(
                        Input(
                            self.base_path,
                            input_path,
                            type,
                            soakdb_path,
                            panddas_paths,
                            excluded_datasets,
                            code_prefix=code_prefix,
                            code_prefix_tooltip=code_prefix_tooltip,
                            logger=self.logger,
                        )
                    )

                elif type == Constants.CONFIG_TYPE_MANUAL:
                    # Determine which datasets to exclude
                    excluded_datasets = utils.find_property(input, Constants.CONFIG_EXCLUDE)
                    if not excluded_datasets:
                        excluded_datasets = []

                    self.logger.info("adding input", input_path)
                    self.inputs.append(
                        Input(
                            self.base_path,
                            input_path,
                            type,
                            None,
                            [],
                            excluded_datasets,
                            code_prefix=code_prefix,
                            code_prefix_tooltip=code_prefix_tooltip,
                            logger=self.logger,
                        )
                    )
                else:
                    raise ValueError("unexpected input type:", type)

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def _log_warning(self, msg):
        self.logger.warn(msg)
        self.warnings.append(msg)

    def _migrate_version(self):
        inp = input("Do you want an environment for a new version of your data creating? (Y/N)")
        if inp == "Y" or inp == "y":
            self.logger.info("migrating data for new data format version", utils.DATA_FORMAT_VERSION)

            new_version_dirname = "upload-v" + str(math.floor(utils.DATA_FORMAT_VERSION))
            new_version_path = self.working_dir / new_version_dirname
            self.logger.info("creating new working dir", new_version_path)
            os.mkdir(new_version_path)
            os.mkdir(new_version_path / "upload_1")
            self.logger.info("copying config.yaml and assemblies.yaml")
            f = shutil.copy2(self.output_path / "config.yaml", new_version_path)
            f = shutil.copy2(self.output_path / "assemblies.yaml", new_version_path)
            extra_files = Path(self.output_path / "extra_files")
            if extra_files.is_dir():
                copy_tree(str(extra_files), str(new_version_path / "extra_files"))
            self.logger.info("removing", self.output_path, "symlink")
            self.output_path.unlink()
            self.logger.info("creating symlink", self.output_path, "->", new_version_path)
            cwd = Path.cwd()
            os.chdir(self.working_dir)
            os.symlink(new_version_dirname, "upload-current", target_is_directory=True)
            os.chdir(cwd)
            self.logger.info(
                "A new directory",
                new_version_dirname,
                "for data format version",
                utils.DATA_FORMAT_VERSION,
                "has been created and the current config.yaml and assemblies.yaml",
                "have been copied there.",
                "\n      It is possible that you might need to update those files.",
                "\n      The old data is in a directory named upload_v? where ? is the old version number",
                "\n      Once ready you can re-run collator using the same command you just used.",
            )
        else:
            self.logger.info(
                "You chose not to create an environment for a new version of the data. "
                + "You will either need to do this manually, or re-run collator and choose to create one"
            )
            exit(0)

    def validate(self):
        v_dir = self.read_versions()
        if not v_dir:
            self.logger.error("Error with version dir. Please fix and try again.")
            return None, None, None
        self.logger.info("using version dir {}".format(v_dir))

        self.logger.info("validating paths")
        self.validate_paths()

        if len(self.logger.errors) == 0:
            self.logger.info("validating data")
            meta = self.validate_data()
        else:
            meta = {}

        return meta

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

    def validate_data(self):
        """
        Read info from the SoakDB database and verify that the necessary entries are present
        :return: The generated metadata
        """

        git_info = None

        if self.include_git_info:
            repo_dir = os.environ.get(Constants.ENV_XCA_GIT_REPO)
            if repo_dir:
                if not Path(repo_dir).is_dir():
                    self._log_error("XCA_GIT_REPO environment variable is defined but the directory does not exist")
            else:
                repo_dir = "./"
            self.logger.info("using GIT repo of", repo_dir)

            try:
                git_info = repo_info.generate(repo_dir)
            except Exception as exc:
                self._log_error(exc)
                self._log_error(traceback.format_exc())
                self._log_error(
                    "cannot determine the status of the Git repo. "
                    + "Is the XCA_GIT_REPO environment variable defined correctly or if not defined is "
                    + "the current directory a Git repo?"
                )

        crystals = {}
        input_dirs = []
        prev_version_dirs_str = [str(d) for d in self.previous_version_dirs]
        cwd = os.getcwd()

        meta = {
            Constants.META_RUN_ON: str(datetime.datetime.now()),
            Constants.CONFIG_CWD: cwd,
            Constants.CONFIG_CONFIG_FILE: str(self.config_file),
            Constants.META_INPUT_DIRS: input_dirs,
            Constants.CONFIG_OUTPUT_DIR: str(self.output_path),
            Constants.META_DATA_FORMAT_VERSION: utils.DATA_FORMAT_VERSION,
            Constants.META_VERSION_NUM: self.version_number,
            Constants.META_VERSION_DIR: str(self.version_dir),
            Constants.META_PREV_VERSION_DIRS: prev_version_dirs_str,
        }
        if git_info:
            meta[Constants.META_GIT_INFO] = git_info

        if self.all_excludes:
            meta[Constants.META_EXClUDES] = list(self.all_excludes)

        tooltips = {}
        for input in self.inputs:
            code_prefix = input.code_prefix
            code_prefix_tooltip = input.code_prefix_tooltip
            if code_prefix and code_prefix_tooltip:
                if code_prefix in tooltips and tooltips[code_prefix] != code_prefix_tooltip:
                    self._log_warning(
                        'code_prefix_tooltip for "{}" is being redefined from "{}" to "{}". To avoid this use unique values for code_prefix.'.format(
                            code_prefix, tooltips[code_prefix], code_prefix_tooltip
                        )
                    )
                tooltips[code_prefix] = code_prefix_tooltip
        if tooltips:
            meta[Constants.META_CODE_PREFIX_TOOLTIPS] = tooltips

        meta[Constants.META_XTALS] = crystals

        for input in self.inputs:
            input_dirs.append(str(input.get_input_dir_path()))
            self._validate_input(input, crystals)

        self._validate_references(crystals)
        self._validate_assemblies(crystals)

        return meta

    def _validate_references(self, crystals):
        refs = utils.find_property(self.config, Constants.CONFIG_REF_DATASETS)
        if refs is None or len(refs) == 0:
            self.logger.info(
                "no references are defined. Use the ref_datasets section of the config if you want to define these"
            )
        else:
            for ref in refs:
                if crystals.get(ref) is None:
                    self._log_error("reference {} is not in the set of crystals to be processed".format(ref))

    def _validate_assemblies(self, crystals):
        assemblies_file = self.output_path / Constants.ASSEMBLIES_FILENAME
        assemblies_yaml = utils.read_config_file(assemblies_file)

        # check the assemblies section
        assemblies = assemblies_yaml.get(Constants.META_ASSEMBLIES)
        if not assemblies:
            self._log_error("assemblies.yaml does not appear to contain any assemblies")
        else:
            for name, data in assemblies.items():
                ref = data[Constants.META_REFERENCE]
                if crystals.get(ref) is None:
                    self._log_error(
                        "reference {} for assembly {} is not in the set of crystals to be processed. Please update your assemblies.yaml file".format(
                            ref, name
                        )
                    )

        # check the crystalforms section
        xtalforms = assemblies_yaml.get(Constants.META_XTALFORMS)
        if not xtalforms:
            self._log_error("assemblies.yaml does not appear to contain any crystalforms")
        else:
            for name, data in xtalforms.items():
                ref = data[Constants.META_REFERENCE]
                if crystals.get(ref) is None:
                    self._log_error(
                        "reference {} for crystalform {} is not in the set of crystals to be processed. Please update your assemblies.yaml file".format(
                            ref, name
                        )
                    )

    def _validate_input(self, input, crystals):
        if input.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            self._validate_soakdb_input(input, crystals)
        elif input.type == Constants.CONFIG_TYPE_MANUAL:
            self._validate_manual_input(input, crystals)
        else:
            raise ValueError("unexpected input type:", input.type)

    def _validate_soakdb_input(self, input, crystals):
        ref_datasets = set(self.config.get(Constants.CONFIG_REF_DATASETS, []))
        dbfile = input.get_soakdb_file_path()
        self.logger.info("opening DB file:", dbfile)
        df = dbreader.filter_dbmeta(dbfile, ref_datasets)
        count = 0
        processed = 0
        num_pdb_files = 0
        num_mtz_files = 0
        num_cif_files = 0

        extra_data = {}

        missing_pdbs = []
        for index, row in df.iterrows():
            count += 1
            xtal_name = row[Constants.SOAKDB_XTAL_NAME]
            cmpd_code = row[Constants.SOAKDB_COL_COMPOUND_CODE]
            cmpd_smiles = row[Constants.SOAKDB_COL_COMPOUND_SMILES]

            # Exclude datasets
            if xtal_name in input.exclude:
                self._log_warning(f"excluding dataset: {xtal_name}")
                continue

            status_str = str(row[Constants.SOAKDB_COL_REFINEMENT_OUTCOME])
            if status_str.startswith("7"):
                # need to check that this crystal was not encountered earlier, if so deprecate it
                self.rejected_xtals.add(xtal_name)
            else:
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
                        inputpath = utils.make_path_relative(Path(file))
                        full_inputpath = self.base_path / inputpath
                        ok = full_inputpath.exists()
                        if ok:
                            num_pdb_files += 1
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            m = (
                                "PDB file for {} not found: {}. "
                                + "Add this to the inputs.exclude section of your config.yaml file if you want to continue"
                            )
                            self._log_error(m.format(xtal_name, full_inputpath))
                            missing_pdbs.append(xtal_name)
                            continue

                        # if we have a PDB file then continue to look for the others
                        colname = Constants.SOAKDB_COL_MTZ
                        file = row[colname]
                        # RefinementMTZ_latest file names are specified as absolute file names, but need to be
                        # handled as relative to the base_path
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
                            self._log_warning(
                                "Crystal {} already exists, it's data will be overriden".format(xtal_name)
                            )

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
                        data[Constants.META_REFINEMENT_OUTCOME] = row[Constants.SOAKDB_COL_REFINEMENT_OUTCOME]
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
                        if cmpd_code:
                            tokens = cmpd_code.split(";")
                            self.compound_codes[xtal_name] = tokens
                        if cmpd_smiles:
                            self.compound_smiles[xtal_name] = utils.parse_compound_smiles(cmpd_smiles)

                        if input.code_prefix is not None:
                            data[Constants.META_CODE_PREFIX] = input.code_prefix
                        data[Constants.META_XTAL_FILES] = f_data

        self.logger.info("validator handled {} rows from database, {} were valid".format(count, processed))
        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding MTZ files".format(num_pdb_files, num_mtz_files)
            )
        if num_cif_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding CIF files".format(num_pdb_files, num_cif_files)
            )
        if missing_pdbs:
            self.logger.warn(
                "PDB files for these crystals were missing. Add them to your inputs.exclude section: "
                + ",".join(missing_pdbs)
            )

    def _validate_manual_input(self, input, crystals):
        num_pdb_files = 0
        num_mtz_files = 0
        num_cif_files = 0
        ref_datasets = set(self.config.get(Constants.CONFIG_REF_DATASETS, []))
        items = self._collect_manual_files(self.base_path / input.input_dir_path)
        self.logger.info("found {} manual inputs".format(len(items)))

        for key, item in items.items():
            if key in input.exclude:
                self.logger.info("excluding manual crystal", key)
                continue
            pdb = item[0]
            mtz = item[1]
            cif = item[2]

            self.logger.info("adding crystal (manual)", pdb.name)
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
            if cif:
                digest = utils.gen_sha256(cif)
                data[Constants.META_XTAL_CIF] = {
                    Constants.META_FILE: cif.relative_to(self.base_path),
                    Constants.META_SHA256: digest,
                }
                num_cif_files += 1
            crystals[key] = {}
            if key in ref_datasets:
                crystals[key][Constants.META_REFERENCE] = True
            if input.code_prefix is not None:
                crystals[key][Constants.META_CODE_PREFIX] = input.code_prefix
            crystals[key][Constants.CONFIG_TYPE] = Constants.CONFIG_TYPE_MANUAL
            crystals[key][Constants.META_XTAL_FILES] = data

        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                "{} PDB files were found, but only {} had corresponding MTZ files".format(num_pdb_files, num_mtz_files)
            )

    def _validate_pdb_file(self, pdb_path, ligand_names):
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("HETATM"):
                    continue

                residue_number = int(line[22:26].strip())
                residue_name = line[17:21].strip()

                if residue_name not in ligand_names:
                    continue

                chain = line[21:22].strip()
                alt_code = line[16:17].strip() or None

                if alt_code:
                    self._log_error(
                        "Encountered ligand "
                        + residue_name
                        + " that has an alternative conformation. "
                        + "You should model this as separate ligands (separate residue numbers)",
                    )

    def _collect_manual_files(self, manual_input_path: Path):
        data = utils.collect_manual_files(manual_input_path)
        self.logger.info(len(data), "manual PDBs found for", manual_input_path)
        return data

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
            self.logger.error(
                "No version directory found. Please create one named", str(self.output_path / "upload_1")
            )
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
            last_meta = self.meta_history[-1]
            last_data_format_version = last_meta.get(utils.Constants.META_DATA_FORMAT_VERSION)
            self.logger.info(
                "Data format versions: current="
                + str(utils.DATA_FORMAT_VERSION)
                + " previous="
                + str(last_data_format_version)
            )
            if last_data_format_version is None:
                self._log_warning(
                    "Previous data format version not defined - cannot determine compatibility"
                    + " with current version"
                )
            elif utils.check_data_format_version(last_data_format_version) < 0:
                self._log_error("Old upload version found that is incompatible with this version of XCA")
                self._migrate_version()
                exit(1)

        return v_dir

    def read_metadata(self, version_dir):
        meta_file = version_dir / Constants.METADATA_XTAL_FILENAME.format("")
        self.logger.info("reading metadata for version {}".format(version_dir), meta_file)

        meta = utils.read_config_file(str(meta_file))
        return meta

    def run(self, meta):
        self.logger.info("running collator...")

        self.logger.info("coping files ...")
        self._write_metadata(meta, suffix="_0")
        new_meta = self._copy_files(meta)
        self.logger.info("munging the history ...")
        self._write_metadata(new_meta, suffix="_1")
        all_xtals, new_xtals = self._munge_history(meta)
        self.logger.info("writing metadata ...")
        self._write_metadata(new_meta)
        self._copy_extra_files()
        self._handle_soakdb()
        self.logger.info("copying config.yaml")
        self._copy_config()
        self.logger.info("run complete")
        return new_meta

    def _handle_soakdb(self):
        names = []
        dataFrames = []
        dbfiles = []
        for input in self.inputs:
            if input.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                dbfile = input.get_soakdb_file_path()
                self.logger.info("reading DB file:", dbfile)
                name, df = dbreader.read_all(dbfile)
                if name in names:
                    self._log_error(
                        "Name {} already exists in a previous input. "
                        + "Make sure the LabVisit column in the soakDB table has distinct values"
                    )
                    return
                names.append(name)
                dataFrames.append(df)
                dbfiles.append(dbfile)
        i = 0
        for name, df, dbfile in zip(names, dataFrames, dbfiles):
            i += 1
            outfile = str(self.output_path / self.version_dir / "extra_files" / ("soakdb_" + name))
            self.logger.info("writing soakdb data to", outfile)
            df.to_csv(outfile + ".csv", index=False)
            f = shutil.copy2(dbfile, outfile + ".sqlite")
            if not f:
                self.logger.error("Failed to copy SoakDB file {}".format(dbfile))

    def _copy_extra_files(self):
        extra_files_dir = self.config.get("extra_files_dir")
        if extra_files_dir is None:
            extra_files_path = self.output_path / "extra_files"
        else:
            extra_files_path = Path(extra_files_dir)

        if not extra_files_path.is_dir():
            self._log_warning("extra_files dir {} not found".format(extra_files_path))
        else:
            src = str(extra_files_path)
            dst = str(self.output_path / self.version_dir / "extra_files")
            self.logger.info("copying extra_files from {} to {}".format(src, dst))
            copied = dir_util.copy_tree(src, dst)
            self.logger.info("copied {} extra_files files".format(len(copied)))

    def _copy_files(self, meta):
        cryst_path = self.version_dir / Constants.META_XTAL_FILES
        ext_cryst_path = self.output_path / cryst_path
        self.logger.info("using cryst_dir of", ext_cryst_path)
        if ext_cryst_path.exists():
            self.logger.info("removing old cryst_dir")
            shutil.rmtree(ext_cryst_path)
        self.logger.info("creating cryst_dir", ext_cryst_path)
        os.makedirs(ext_cryst_path)

        extra_files_path = self.output_path / self.version_dir / "extra_files"
        if not extra_files_path.is_dir():
            os.mkdir(extra_files_path)

        compound_path = extra_files_path / "compounds_auto.csv"
        self.logger.info("writing compound data to " + str(compound_path))
        with open(compound_path, "wt") as compounds_auto:
            # write header for compounds_auto.csv
            compounds_auto.write("xtal,ligand_name," + Constants.META_CMPD_CODE + "\n")

            event_tables = self._find_event_tables()
            forbidden_unattested_ligand_events = {}
            pdb_count = 0
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
                            # PDB is new or changed
                            if old_digest is not None:
                                self.logger.info("PDB file for {} has changed".format(xtal_name))
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
                            old_digest = historical_xtal_data.get(Constants.META_XTAL_MTZ, {}).get(
                                Constants.META_SHA256
                            )
                            if digest != old_digest:
                                if old_digest is not None:
                                    self.logger.info("MTZ file for {} has changed".format(xtal_name))
                                mtz_name = xtal_name + ".mtz"
                                mtz_output = dir / mtz_name
                                files_to_copy[Constants.META_XTAL_MTZ] = (mtz_input, mtz_output, digest)
                        elif type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                            self.logger.warn("mtz file {} not present".format(mtz_input))
                    else:
                        self.logger.warn("MTZ entry missing for {}".format(xtal_name))

                    # handle the CIF file
                    cif = curr_xtal_data.get(Constants.META_XTAL_CIF)
                    ligand_names = []
                    ligand_mols = []
                    if cif:
                        cif_file = cif[Constants.META_FILE]
                        cif_input = self.base_path / cif_file
                        if cif_input.is_file():
                            digest = utils.gen_sha256(cif_input)
                            old_digest = historical_xtal_data.get(Constants.META_XTAL_CIF, {}).get(
                                Constants.META_SHA256
                            )
                            if digest != old_digest:
                                if old_digest is not None:
                                    self.logger.info("CIF file for {} has changed".format(xtal_name))
                                cif_name = xtal_name + ".cif"
                                cif_output = dir / cif_name
                                files_to_copy[Constants.META_XTAL_CIF] = (cif_input, cif_output, digest)
                            # read the CIF file to determine the ligand names
                            try:
                                ligand_mols = utils.gen_mols_from_cif(cif_input)
                                ligand_names = [m.GetProp("_Name") for m in ligand_mols]
                                self.logger.info("Found ligand names " + str(ligand_names))
                            except:
                                tb = traceback.format_exc()
                                self._log_error(
                                    "Failed to handle CIF file " + str(cif_input) + " for crystal " + xtal_name
                                )
                                self.logger.info("Error was:\n" + tb)
                                continue

                    elif type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                        self.logger.warn("CIF entry missing for {}".format(xtal_name))

                    # Handle historical ligand binding events (in particular pull up their event map SHA256s for comparing)
                    hist_event_maps = {}
                    for ligand_binding_data in historical_xtal_data.get(Constants.META_BINDING_EVENT, []):
                        model = ligand_binding_data.get(Constants.META_PROT_MODEL)
                        chain = ligand_binding_data.get(Constants.META_PROT_CHAIN)
                        res = ligand_binding_data.get(Constants.META_PROT_RES)
                        if model is not None and chain and res:
                            hist_event_maps[(model, chain, res)] = ligand_binding_data

                    # Determine the ligands present and their coordinates
                    dataset_ligands = self.get_dataset_ligands(pdb_input, ligand_names)
                    if xtal_name == "Mpro-i0130":
                        print(f"pdb_input: {pdb_input}")
                        print(f"Ligand names: {ligand_names}")
                        print(f"Got ligands: {dataset_ligands}")

                    # Match ligand to panddas event maps if possible and determine if those maps are new
                    best_event_map_paths = self.get_dataset_event_maps(xtal_name, dataset_ligands, event_tables)
                    identical_historical_event_maps = {}
                    unattested_ligand_events = {}
                    attested_ligand_events = {}
                    event_maps_to_copy = {}

                    for ligand_key in dataset_ligands:
                        if ligand_key in best_event_map_paths:
                            ligand_event_map_data = best_event_map_paths[ligand_key]
                            path = ligand_event_map_data[0]
                            if path:
                                digest = utils.gen_sha256(path)
                                ccp4_output = (
                                    cryst_path
                                    / xtal_name
                                    / "{}_{}_{}_{}.ccp4".format(xtal_name, ligand_key[0], ligand_key[1], ligand_key[2])
                                )
                                attested_ligand_events[ligand_key] = (
                                    path,
                                    ccp4_output,
                                    digest,
                                    ligand_key,
                                    ligand_event_map_data[1],
                                    ligand_event_map_data[2],
                                )
                                hist_data = hist_event_maps.get(ligand_key)
                                # Track whether the event map actually is new by data
                                if hist_data:
                                    if digest == hist_data.get(Constants.META_SHA256):
                                        identical_historical_event_maps[ligand_key] = True
                                        self.logger.info(
                                            "Event map {} for {} is unchanged".format(ligand_key, xtal_name)
                                        )
                                    else:
                                        event_maps_to_copy[ligand_key] = True
                                        self.logger.info(
                                            "Event map {} for {} has changed".format(ligand_key, xtal_name)
                                        )
                                else:
                                    event_maps_to_copy[ligand_key] = True
                        # Handle ligands that cannot be matched
                        else:
                            # Add those permitted ligands
                            if xtal_name in self.panddas_missing_ok:
                                self.logger.warn(
                                    "no PanDDA event map found for",
                                    xtal_name,
                                    "but this is OK as it's been added to the panddas_missing_ok",
                                    "list in the config file",
                                )
                                unattested_ligand_events[ligand_key] = True
                            # Track forbidden ligands for informative error messages at the end of this function
                            else:
                                self.logger.error(
                                    "no PanDDA event map found. If you want to allow this then add",
                                    xtal_name,
                                    "to the panddas_missing_ok list in the config file",
                                )
                                forbidden_unattested_ligand_events[xtal_name] = ligand_key

                else:
                    self.logger.error("PDB entry missing for {}".format(xtal_name))
                    return meta

                # now copy the files
                self.logger.info("{} has {} files to copy".format(xtal_name, len(files_to_copy)))
                fdata = files_to_copy.get(Constants.META_XTAL_PDB)
                data_to_add = {}
                if fdata:
                    os.makedirs(self.output_path / dir)
                    f = shutil.copy2(fdata[0], self.output_path / fdata[1])
                    if not f:
                        self.logger.error(
                            "Failed to copy PDB file {} to {}".format(fdata[0], self.output_path / fdata[1])
                        )
                    else:
                        pdb_count += 1
                        data_to_add[Constants.META_XTAL_PDB] = {
                            Constants.META_FILE: str(fdata[1]),
                            Constants.META_SHA256: fdata[2],
                            Constants.META_SOURCE_FILE: str(fdata[0]),
                        }
                        # copy MTZ file
                        fdata = files_to_copy.get(Constants.META_XTAL_MTZ)
                        if fdata:
                            f = shutil.copy2(fdata[0], self.output_path / fdata[1])
                            if not f:
                                self.logger.error(
                                    "Failed to copy MTZ file {} to {}".format(fdata[0], self.output_path / fdata[1])
                                )
                            else:
                                data_to_add[Constants.META_XTAL_MTZ] = {
                                    Constants.META_FILE: str(fdata[1]),
                                    Constants.META_SHA256: fdata[2],
                                    Constants.META_SOURCE_FILE: str(fdata[0]),
                                }
                        fdata = files_to_copy.get(Constants.META_XTAL_CIF)

                        # copy CIF file
                        if fdata:
                            f = shutil.copy2(fdata[0], self.output_path / fdata[1])
                            if not f:
                                self.logger.error(
                                    "Failed to copy CIF file {} to {}".format(fdata[0], self.output_path / fdata[1])
                                )
                            else:
                                data_to_add[Constants.META_XTAL_CIF] = {
                                    Constants.META_FILE: str(fdata[1]),
                                    Constants.META_SHA256: fdata[2],
                                    Constants.META_SOURCE_FILE: str(fdata[0]),
                                }
                                try:
                                    ligands = {}

                                    cpd_codes = self.compound_codes.get(xtal_name)
                                    if cpd_codes and len(cpd_codes) == len(ligand_mols):
                                        cpd_codes_is_valid = True
                                    else:
                                        cpd_codes_is_valid = False
                                        if cpd_codes:
                                            self._log_error("Invalid number of compound codes for " + xtal_name)
                                        # else:  no compound codes defined - this is OK

                                    cpd_smiles = self.compound_smiles.get(xtal_name)
                                    if cpd_smiles and len(cpd_smiles) == len(ligand_mols):
                                        cpd_smiles_is_valid = True
                                    else:
                                        cpd_smiles_is_valid = False
                                        if cpd_smiles:
                                            self._log_error("Invalid number of compound SMILES for " + xtal_name)
                                        # else:  no compound smiles defined - this is OK

                                    for i, mol in enumerate(ligand_mols):
                                        name = mol.GetProp("_Name")
                                        smi = None
                                        ligands[name] = {}
                                        if cpd_codes_is_valid:
                                            ligands[name][Constants.META_CMPD_CODE] = cpd_codes[i]
                                            compounds_auto.write(",".join((xtal_name, name, cpd_codes[i])) + "\n")
                                        if cpd_smiles_is_valid:
                                            ligands[name][Constants.META_MODELED_SMILES_SOAKDB] = cpd_smiles[i][0]
                                            try:
                                                m = Chem.MolFromSmiles(cpd_smiles[i][0])
                                                can_smi = Chem.MolToSmiles(m)
                                                if can_smi:
                                                    ligands[name][Constants.META_MODELED_SMILES_CANON] = can_smi
                                                    smi = can_smi
                                                else:
                                                    self._log_warning(
                                                        "Could not generate canonical SMILES for the modelled SMILES "
                                                        + "in SoakDB - instead using the non-canonical form"
                                                    )
                                                    smi = cpd_smiles[i][0]
                                            except:
                                                self._log_warning(
                                                    "Failed to generate canonical smiles for "
                                                    + "modeled molecule from soakDB"
                                                )
                                            if len(cpd_smiles[i]) == 2:
                                                ligands[name][Constants.META_SOAKED_SMILES_SOAKDB] = cpd_smiles[i][1]
                                                try:
                                                    m = Chem.MolFromSmiles(cpd_smiles[i][1])
                                                    can_smi = Chem.MolToSmiles(m)
                                                    if can_smi:
                                                        ligands[name][Constants.META_MODELED_SMILES_CANON] = can_smi
                                                        smi = can_smi
                                                    else:
                                                        smi = cpd_smiles[i][1]
                                                except:
                                                    self._log_warning(
                                                        "Failed to generate canonical smiles for "
                                                        + "soaked molecule from soakDB"
                                                    )
                                        if not smi:
                                            try:
                                                smi = Chem.MolToSmiles(mol)
                                            except:
                                                self._log_warning(
                                                    "Failed to generate SMILES from CIF molecule for ligand " + name
                                                )
                                        if smi:
                                            ligands[name][Constants.META_SMILES] = smi
                                        else:
                                            self._log_error(
                                                "could not generate SMILES for "
                                                + xtal_name
                                                + " - not defined in SoakDB CompoundSMILES column "
                                                + "and could not be generated from CIF file"
                                            )

                                    if ligands:
                                        data_to_add[Constants.META_XTAL_CIF][Constants.META_LIGANDS] = ligands
                                except:
                                    self.logger.warn("failed to generate ligand data for {}".format(xtal_name))
                                    traceback.print_exc()

                        # copy event maps that differ in SHA from previously known ones
                        unsucessfully_copied_event_maps = {}
                        if len(event_maps_to_copy) != 0:
                            for ligand_key in event_maps_to_copy:
                                source = attested_ligand_events[ligand_key][0]
                                destination = attested_ligand_events[ligand_key][1]
                                f = shutil.copy2(source, self.output_path / destination)
                                if not f:
                                    self.logger.error(
                                        "Failed to copy Panddas file {} to {}".format(
                                            source, self.output_path / destination
                                        )
                                    )
                                    # Mark that copying failed
                                    unsucessfully_copied_event_maps[ligand_key] = True

                        # Create ligand binding events for the dataset
                        ligand_binding_events = []
                        for ligand_key in dataset_ligands:
                            # Add binding events for ligands that can be matched to PanDDA event maps
                            if ligand_key in attested_ligand_events:
                                # Skip if failed to copy pandda event map
                                if ligand_key in unsucessfully_copied_event_maps:
                                    continue
                                attested_ligand_event_data = attested_ligand_events[ligand_key]

                                if identical_historical_event_maps.get(ligand_key, False):
                                    data = hist_event_maps[ligand_key]
                                else:
                                    data = {
                                        Constants.META_FILE: str(attested_ligand_event_data[1]),
                                        Constants.META_SHA256: attested_ligand_event_data[2],
                                        Constants.META_SOURCE_FILE: str(attested_ligand_event_data[0]),
                                        Constants.META_PROT_MODEL: ligand_key[0],
                                        Constants.META_PROT_CHAIN: ligand_key[1],
                                        Constants.META_PROT_RES: ligand_key[2],
                                        Constants.META_PROT_NAME: ligand_key[3],
                                        Constants.META_PROT_INDEX: attested_ligand_event_data[4],
                                        Constants.META_PROT_BDC: attested_ligand_event_data[5],
                                    }
                                # if data[Constants.META_FILE]:
                                ligand_binding_events.append(data)
                            # Add binding events for permitted ligands without an event map
                            elif ligand_key in unattested_ligand_events:
                                data = {
                                    Constants.META_PROT_MODEL: ligand_key[0],
                                    Constants.META_PROT_CHAIN: ligand_key[1],
                                    Constants.META_PROT_RES: ligand_key[2],
                                }
                                ligand_binding_events.append(data)

                            # Skip if ligand key is not associated with a legal ligand
                            else:
                                continue
                            # ligand_binding_events.append(data)

                        # Add data on the ligand binding events to the new dataset to add
                        if ligand_binding_events:
                            data_to_add[Constants.META_BINDING_EVENT] = ligand_binding_events

                self.num_crystals = pdb_count

                new_xtal_data = {}
                for k, v in historical_xtal_data.items():
                    new_xtal_data[k] = v
                for k, v in data_to_add.items():
                    new_xtal_data[k] = v
                xtal[Constants.META_XTAL_FILES] = new_xtal_data

            # Handle the presence of ligand without event maps that have not been permitted
            if len(forbidden_unattested_ligand_events) != 0:
                exception = (
                    "No PanDDA event map found that correspond to the following ligands. If you want to allow these then "
                    "add the corresponding crystal names to the panddas_missing_ok list in the config file:\n"
                )
                for dtag, ligand_key in forbidden_unattested_ligand_events.items():
                    lk = ligand_key
                    exception = exception + f"{dtag} : Model: {lk[0]}; Chain: {lk[1]}; Residue: {lk[2]}\n"
                raise Exception(exception)

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
        overrides = self.config.get(Constants.CONFIG_OVERRIDES, {})

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
        # if status is now '7 - Analysed & Rejected' then it will be in the set of rejected_xtals so
        # we must mark the status as deprecated
        if xtal_name in self.rejected_xtals:
            return Constants.META_STATUS_DEPRECATED

        old_pdb_data = old_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        new_pdb_data = new_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        old_sha256 = old_pdb_data[Constants.META_SHA256]
        new_sha256 = new_pdb_data[Constants.META_SHA256]
        if old_sha256 == new_sha256:
            return Constants.META_STATUS_UNCHANGED
        else:
            return Constants.META_STATUS_SUPERSEDED

    def _write_metadata(self, meta, suffix=""):
        f = self.output_path / self.version_dir / Constants.METADATA_XTAL_FILENAME.format(suffix)
        with open(f, "w") as stream:
            yaml.dump(meta, stream, sort_keys=False)

    def _copy_config(self):
        to_path = self.output_path / self.version_dir / "config.yaml"
        self.logger.info("copying config file", self.config_file, "to", to_path)
        f = shutil.copy2(self.config_file, to_path)
        if not f:
            self.logger.warn("Failed to copy config file to {}".format((self.output_path / self.version_dir)))
            return False
        return True

    def get_dataset_ligands(
        self, structure_path: Path, ligand_names: list[str]
    ) -> dict[tuple[str, str, str], np.array]:
        structure = gemmi.read_structure(str(structure_path))

        ligand_coords = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name in ligand_names:
                        poss = []
                        for atom in residue:
                            pos = atom.pos
                            poss.append([pos.x, pos.y, pos.z])

                        arr = np.array(poss)
                        mean = np.mean(arr, axis=0)
                        ligand_coords[(model.name, chain.name, residue.seqid.num, residue.name)] = mean

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
        self, xtal_name: str, ligand_coords: dict[(str, str, int), np.array], event_tables: dict[Path, pd.DataFrame]
    ) -> dict[tuple[str, str, str], Path]:
        # Get the relevant structure

        # Get the coordinates of ligands
        # ligand_coords = self.get_ligand_coords(structure)

        # Get the closest events within some reasonable radius
        closest_event_maps = {}
        for ligand_key, ligand_coord in ligand_coords.items():
            closest_event_map, data = self.get_closest_event_map(xtal_name, ligand_coord, event_tables)
            if closest_event_map:
                closest_event_maps[ligand_key] = (closest_event_map, data[0], data[1])

        return closest_event_maps


def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-d", "--dir", help="Working directory")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("-v", "--validate", action="store_true", help="Only perform validation")
    parser.add_argument("--no-git-info", action="store_false", help="Don't add GIT info to metadata")

    args = parser.parse_args()

    if args.dir:
        working_dir = Path(args.dir)
    else:
        working_dir = Path.cwd()

    wd = utils._verify_working_dir(working_dir)

    if wd:
        working_dir = wd
    else:
        print("Working dir does not seem to have been initialised - missing 'upload_current' symlink")
        inp = input("Do you want the working dir to be initialised? (Y/N)")
        if inp == "Y" or inp == "y":
            print("Initialising working dir")
            s = setup.Setup(args.dir)
            s.run()
        exit(1)

    logger = None
    try:
        c = Collator(working_dir, log_level=args.log_level, include_git_info=args.no_git_info)

        logger = c.logger
        logger.info("collator: ", str(args))
        utils.LOG = logger

        meta = c.validate()

        if not args.validate:
            if meta is None or len(logger.errors) > 0:
                logger.error("There are errors, cannot continue")
                logger.report()
                logger.close()
                exit(1)
            else:
                t0 = time.time()
                c.run(meta)
                t1 = time.time()
                # write a summary of errors and warnings
                logger.info("Handled {} crystals in {} secs".format(c.num_crystals, round(t1 - t0)))
                logger.report()
                logger.close()
                if logger.logfilename:
                    to_path = c.output_path / c.version_dir / "collator.log"
                    print("copying log file", logger.logfilename, "to", to_path)
                    f = shutil.copy2(logger.logfilename, to_path)
                    if not f:
                        print("Failed to copy log file {} to {}".format(logger.logfilename, to_path))

    except Exception as err:
        # uncaught exception
        tb = traceback.format_exc()
        if logger:
            logger.error(
                "Unexpected fatal error occurred when running collator\n"
                + tb
                + "\nPlease send this information to the XCA developers:\nLog file: "
                + str(logger.logfilename)
                + "\nWorking dir location: "
                + str(wd)
            )
        else:  # couldn't even create the Collator object so log file will contain nothing useful
            print(tb)
            print(
                "Unexpected fatal error occurred, most likely the configuration is wrong or collator was invoked "
                + "incorrectly. Please send the command you ran and your current directory to the developers."
            )


if __name__ == "__main__":
    main()
