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

import datetime
import re
from pathlib import Path

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


def expand_file_path(filepath: Path, default='.'):
    if filepath.is_absolute():
        return filepath
    else:
        return Path(default) / filepath


class Input:

    def __init__(self, base_path: Path, input_dir_path: Path, type: str, soakdb_file_path,
                 panddas_event_file_paths: list[Path]):
        self.base_path = base_path
        self.input_dir_path = input_dir_path
        self.type = type
        self.soakdb_file_path = soakdb_file_path
        self.panddas_event_file_paths = panddas_event_file_paths
        self.errors = []
        self.warnings = []

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
            self.errors.append('base path must be defined')
        elif not self.base_path.exists():
            self.errors.append('base_path does not exist: {}'.format(self.base_path))
        elif not self.base_path.is_dir():
            self.errors.append('base_path is not a directory: {}'.format(self.base_path))

        if not self.input_dir_path:
            self.errors.append('input_dir_path is not defined')
        elif self.input_dir_path.is_absolute():
            self.errors.append('input_path must be a path relative to base_path')
        else:
            p = self.get_input_dir_path()
            print('testing', p)
            if not p.exists():
                self.errors.append('input_dir_path does not exist: {}'.format(p))
            elif not p.is_dir():
                self.errors.append('input_dir_path is not a directory: {}'.format(p))

        if self.type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            if not self.soakdb_file_path:
                self.errors.append('soakdb_file_path is not defined')
            elif self.soakdb_file_path.is_absolute():
                self.errors.append('soakdb_file_path must be a path relative to input_path')
            else:
                p = self.get_soakdb_file_path()
                if not p.exists():
                    self.errors.append('soakdb_file_path does not exist: {}'.format(p))
                elif not p.is_file():
                    self.errors.append('soakdb_file_path is not a file: {}'.format(p))

        return len(self.errors), len(self.warnings)


class Processor:
    """
    Base class that understands versioning and basic configuration
    """

    def __init__(self, config_file: Path, logger=None):
        self.errors = []
        self.warnings = []
        self.config_file = config_file

        config = utils.read_config_file(config_file)

        self.base_path = utils.find_path(config, Constants.CONFIG_BASE_DIR)
        self.output_path = utils.find_path(config, Constants.CONFIG_OUTPUT_DIR)

        self.target_name = utils.find_property(config, Constants.CONFIG_TARGET_NAME)
        self.config = config
        self.version_dir = None
        self.meta_history = []
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

        self.inputs = []
        inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
        self.logger.info('found {} inputs'.format(len(inputs)))
        if inputs:
            for input in inputs:
                input_path = utils.find_path(input, Constants.CONFIG_DIR)
                type = utils.find_property(input, Constants.CONFIG_TYPE)
                if type == Constants.CONFIG_TYPE_MODEL_BUILDING:
                    soakdb_path = utils.find_path(input, Constants.CONFIG_SOAKDB, default=Constants.DEFAULT_SOAKDB_PATH)
                    panddas_csvs = utils.find_property(input, Constants.META_BINDING_EVENT)
                    if panddas_csvs:
                        panddas_paths = [Path(p) for p in panddas_csvs]
                    else:
                        panddas_paths = []

                    self.logger.info('adding input', input_path)
                    self.inputs.append(Input(self.base_path, input_path, type, soakdb_path, panddas_paths))

                elif type == Constants.CONFIG_TYPE_MANUAL:
                    self.logger.info('adding input', input_path)
                    self.inputs.append(Input(self.base_path, input_path, type, None, []))

                else:
                    raise ValueError('unexpected input type:', type)

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def _log_warning(self, msg):
        self.logger.warn(msg)
        self.warnings.append(msg)

    def validate_paths(self):
        self.errors.clear()
        self.warnings.clear()

        if not self.target_name:
            self._log_error('target_name not defined')
        else:
            if len(self.target_name) < 4:
                self._log_error('target_name must have at least 4 characters: ' + self.target_name)
            else:
                x = re.search("^[A-Za-z]+[A-Za-z0-9_\\-]*$", self.target_name)
                if not x:
                    self._log_error('Invalid target_name: ' + self.target_name)

        if self.output_path:
            if not self.output_path.exists():
                self._log_error('output_dir does not exist: {}'.format(self.output_path))
            elif not self.output_path.is_dir():
                self._log_error('output_dir argument is not a directory: {}'.format(self.output_path))
        else:
            self._log_error('Output dir not defined in the config file')

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
        meta = {
            Constants.META_RUN_ON: str(datetime.datetime.now()),
            Constants.META_INPUT_DIRS: input_dirs,
            Constants.CONFIG_OUTPUT_DIR: str(self.output_path),
            Constants.META_XTALS: crystals}

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
            raise ValueError('unexpected input type:', input.type)

    def _validate_soakdb_input(self, input, crystals):
        dbfile = input.get_soakdb_file_path()
        self.logger.info('Opening DB file:', dbfile)
        df = dbreader.filter_dbmeta(dbfile)
        count = 0
        processed = 0
        num_pdb_files = 0
        num_mtz_files = 0
        num_cif_files = 0

        for index, row in df.iterrows():
            count += 1
            xtal_name = row[Constants.SOAKDB_XTAL_NAME]
            xtal_dir = generate_xtal_dir(input.input_dir_path, xtal_name)
            # print('xtal_dir:', xtal_dir)
            if not xtal_name:
                self._log_error('Crystal name not defined, cannot process row {}'.format(xtal_name))
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
                    self._log_warning('PDB entry {} for {} not defined in SoakDB'.format(colname, xtal_name))
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
                        self._log_warning('PDB file for {} not found: {}'.format(xtal_name, full_inputpath))

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
                            self._log_warning('MTZ file for {} not found: {}'.format(xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning('MTZ entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    colname = Constants.SOAKDB_COL_CIF
                    file = row[colname]
                    # RefinementCIF file names are relative to the xtal_dir
                    if file:
                        p = Path(file)
                        if p.is_absolute():
                            inputpath = p.relative_to('/')
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
                            self._log_warning('CIF file for {} not found: {}'.format(xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning('CIF entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    processed += 1
                    if xtal_name in crystals.keys():
                        self._log_warning(
                            "Crystal {} already exists, it's data will be overriden".format(xtal_name))

                    data = {Constants.CONFIG_TYPE: Constants.CONFIG_TYPE_MODEL_BUILDING}
                    self.logger.info('adding crystal (model_building)', xtal_name)
                    crystals[xtal_name] = data
                    last_updated_date = row[Constants.SOAKDB_COL_LAST_UPDATED]
                    if last_updated_date:
                        dt_str = last_updated_date.strftime(utils._DATETIME_FORMAT)
                        data[Constants.META_LAST_UPDATED] = dt_str
                    digest = utils.gen_sha256(self.base_path / expanded_files[0])
                    f_data = {Constants.META_XTAL_PDB: {
                        Constants.META_FILE: str(expanded_files[0]),
                        Constants.META_SHA256: digest}}
                    if expanded_files[1]:
                        digest = utils.gen_sha256(self.base_path / expanded_files[1])
                        f_data[Constants.META_XTAL_MTZ] = {
                            Constants.META_FILE: str(expanded_files[1]),
                            Constants.META_SHA256: digest
                        }
                    if expanded_files[2]:
                        digest = utils.gen_sha256(self.base_path / expanded_files[2])
                        f_data[Constants.META_XTAL_CIF] = {
                            Constants.META_FILE: str(expanded_files[2]),
                            Constants.META_SHA256: digest
                        }
                    data[Constants.META_XTAL_FILES] = f_data

        self.logger.info('Validator handled {} rows from database, {} were valid'.format(count, processed))
        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                '{} PDB files were found, but only {} had corresponding MTZ files'.format(num_pdb_files, num_mtz_files))
        if num_cif_files < num_pdb_files:
            self.logger.warn(
                '{} PDB files were found, but only {} had corresponding CIF files'.format(num_pdb_files, num_cif_files))

    def _validate_manual_input(self, input, crystals):
        num_pdb_files = 0
        num_mtz_files = 0
        for child in (self.base_path / input.input_dir_path).iterdir():
            pdb = None
            mtz = None
            if child.is_dir():
                for file in child.iterdir():
                    if file.suffix == '.pdb':
                        pdb = file
                    if file.suffix == '.mtz':
                        mtz = file
                if pdb:
                    self.logger.info('adding crystal (manual)', child.name)
                    num_pdb_files += 1
                    digest = utils.gen_sha256(pdb)
                    data = {Constants.META_XTAL_PDB: {
                        Constants.META_FILE: pdb.relative_to(self.base_path),
                        Constants.META_SHA256: digest}}
                    if mtz:
                        digest = utils.gen_sha256(mtz)
                        data[Constants.META_XTAL_MTZ] = {
                            Constants.META_FILE: mtz.relative_to(self.base_path),
                            Constants.META_SHA256: digest}
                        num_mtz_files += 1
                    crystals[child.name] = {
                        Constants.CONFIG_TYPE: Constants.CONFIG_TYPE_MANUAL,
                        Constants.META_XTAL_FILES: data}

        if num_mtz_files < num_pdb_files:
            self.logger.warn(
                '{} PDB files were found, but only {} had corresponding MTZ files'.format(num_pdb_files, num_mtz_files))

    def read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = self.output_path / Constants.VERSION_DIR_PREFIX + str(version)
            if v_dir.is_dir():
                version += 1
            else:
                break
        if version == 1:
            self.logger.error('No version directory found. Please create one named upload_1')
            return None

        # the working version dir is one less than the current value
        version -= 1
        self.logger.info('Version is {}'.format(version))
        v_dir = self.output_path / Constants.VERSION_DIR_PREFIX + str(version)

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(1, version):
                self.logger.info('Reading metadata for version {}'.format(v))
                dir_name = self.output_path / Constants.VERSION_DIR_PREFIX + str(v)
                meta = self.read_metadata(dir_name)
                self.meta_history.append(meta)

        self.logger.info('Setting version dir to {}'.format(v_dir))
        self.version_dir = v_dir

        num_old_metas = len(self.meta_history)
        if num_old_metas:
            self.logger.info('Found {} metadata files from previous versions'.format(num_old_metas))

        return v_dir

    def read_metadata(self, version_dir):
        self.logger.info('Reading metadata for version {}'.format(version_dir))
        meta_file = version_dir / Constants.METADATA_XTAL_FILENAME

        meta = utils.read_config_file(meta_file)
        return meta
