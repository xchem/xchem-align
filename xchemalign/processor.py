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


from pathlib import Path
import re, datetime

from . import utils, dbreader


METADATA_FILENAME = 'metadata.yaml'
CONFIG_FILENAME = 'config.yaml'
VERSION_DIR_PREFIX = 'upload_'

def generate_xtal_dir(input_path: Path, xtal_name: str):
    """
    Generate the directory with the crystal data
    :param input_path:
    :param xtal_name:
    :return: The Path to the base dir for this crystal
    """
    xtal_dir = input_path / 'processing/analysis/model_building' / xtal_name
    return xtal_dir

def generate_file_paths(filepath: Path, xtal_dir_path: Path, output_path: Path):

    if filepath.is_absolute():
        inputpath = filepath
        outputpath = output_path / filepath.relative_to('/')
    else:
        inputpath = xtal_dir_path / filepath
        outputpath = output_path / xtal_dir_path / filepath

    return inputpath, outputpath

def expand_file_path(filepath: Path, default='.'):
    if filepath.is_absolute():
        return filepath
    else:
        return Path(default) / filepath


class Input:

    def __init__(self, base_path: Path, input_dir_path: Path, soakdb_file_path, panddas_event_file_paths: list[Path]):
        self.base_path = base_path
        self.input_dir_path = input_dir_path
        self.soakdb_file_path = soakdb_file_path
        self.panddas_event_file_paths = panddas_event_file_paths
        self.errors = []
        self.warnings = []

    def expand_path(self, p, expand=True):
        return utils.expand_path(self.base_path, p, expand=expand)


    def get_input_dir_path(self, expand=True):
        return self.expand_path(self.input_dir_path, expand=expand)

    def get_soakdb_file_path(self, expand=True):
        return self.expand_path(self.soakdb_file_path, expand=expand)

    def validate(self):
        if not self.expand_path(self.input_dir_path):
            self.errors.append('input_dir_path is not defined')
        elif not self.expand_path(self.input_dir_path).exists():
            self.errors.append('input_dir_path does not exist: {}'.format(self.input_dir_path))
        elif not self.expand_path(self.input_dir_path).is_dir():
            self.errors.append('input_dir_path is not a directory: {}'.format(self.input_dir_path))

        if not self.expand_path(self.soakdb_file_path):
            self.errors.append('soakdb_file_path is not defined')
        elif not self.expand_path(self.soakdb_file_path).exists():
            self.errors.append('soakdb_file_path does not exist: {}'.format(self.soakdb_file_path))
        elif not self.expand_path(self.soakdb_file_path).is_file():
            self.errors.append('soakdb_file_path is not a file: {}'.format(self.soakdb_file_path))

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

        self.base_path = utils.find_path(config, 'base_dir')
        self.output_path = utils.find_path(config, 'output_dir')
        self.copy_path = utils.find_path(config, 'copy_dir')

        self.target_name =  utils.find_property(config, 'target_name')
        self.config = config
        self.version_dir = None
        self.meta_history = []
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

        self.inputs = []
        inputs = utils.find_property(config, 'inputs')
        if inputs:
            for input in inputs:
                input_path = utils.find_path(input, 'dir')
                soakdb_path = utils.find_path(input, 'soakdb', default='processing/database/soakDBDataFile.sqlite')

            # TODO - handle panddas
            panddas_paths = []

            self.inputs.append(Input(self.base_path, input_path, input_path / soakdb_path, panddas_paths))

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

    def validate_soakdb_data(self):
        """
        Read info from the SoakDB database and verify that the necessary entries are present
        :return: The generated metadata
        """
        valid_ids = {}
        meta = {
            'run_on': str(datetime.datetime.now()),
            'input_dirs': '',
            'output_dir': str(self.output_path),
            'crystals': valid_ids}

        for input in self.inputs:
            dbfile = input.get_soakdb_file_path()
            self.logger.info('Opening DB file:', dbfile)
            df = dbreader.filter_dbmeta(dbfile)
            count = 0
            processed = 0
            for index, row in df.iterrows():
                count += 1
                xtal_name = row['CrystalName']
                xtal_dir = generate_xtal_dir(input.input_dir_path, xtal_name)
                # print('xtal_dir:', xtal_dir)
                if not xtal_name:
                    self._log_error('Crystal name not defined, cannot process row {}'.format(xtal_name))
                else:
                    self.logger.info('Processing crystal {} {}'.format(count, xtal_name))
                    missing_files = 0
                    expanded_files = []
                    colname = 'RefinementPDB_latest'
                    file = row[colname]
                    if file:
                        # print('handling', colname, file)
                        inputpath = expand_file_path(Path(file))
                        full_inputpath = utils.expand_path(self.base_path, inputpath)
                        # print('generated', full_inputpath)
                        ok = full_inputpath.exists()
                        if ok:
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            self._log_warning('File {} for {} not found: {}'.format(colname, xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning('Entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    colname = 'RefinementMTZ_latest'
                    file = row[colname]
                    if file:
                        # print('handling', colname, file)
                        inputpath = expand_file_path(Path(file))
                        full_inputpath = utils.expand_path(self.base_path, inputpath)
                        # print('generated', full_inputpath)
                        ok = full_inputpath.exists()
                        if ok:
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            self._log_warning('File {} for {} not found: {}'.format(colname, xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning('Entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    colname = 'RefinementCIF'
                    file = row[colname]
                    if file:
                        # print('handling', colname, file)
                        inputpath = expand_file_path(Path(file), default=xtal_dir)
                        full_inputpath = utils.expand_path(self.base_path, inputpath)
                        # print('generated', full_inputpath)
                        ok = full_inputpath.exists()
                        if ok:
                            expanded_files.append(inputpath)
                        else:
                            expanded_files.append(None)
                            missing_files += 1
                            self._log_warning('File {} for {} not found: {}'.format(colname, xtal_name, full_inputpath))
                    else:
                        expanded_files.append(None)
                        self._log_warning('Entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    if missing_files > 0:
                        self._log_warning('{} files for {} missing. Will not process'.format(missing_files, xtal_name))
                    else:
                        processed += 1
                        if xtal_name in valid_ids.keys():
                            self._log_warning("Crystal {} already exists, it's data will be overriden".format(xtal_name))

                        data = {}
                        valid_ids[xtal_name] = data
                        last_updated_date = row['LastUpdatedDate']
                        if last_updated_date:
                            dt_str = last_updated_date.strftime(utils._DATETIME_FORMAT)
                            print('date', dt_str)
                            data['last_updated'] = dt_str
                        data['crystallographic_files'] = {
                            'xtal_pdb': expanded_files[0],
                            'xtal_mtz': expanded_files[1],
                            'ligand_cif': expanded_files[2]}

            self.logger.info('Validator handled {} rows from database, {} were valid'.format(count, processed))

        return meta

    def read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = self.output_path / VERSION_DIR_PREFIX + str(version)
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
        v_dir = self.output_path / VERSION_DIR_PREFIX + str(version)

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(1, version):
                self.logger.info('Reading metadata for version {}'.format(v))
                dir_name = self.output_path / VERSION_DIR_PREFIX + str(v)
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
        meta_file = version_dir / METADATA_FILENAME

        meta = utils.read_config_file(meta_file)
        return meta
