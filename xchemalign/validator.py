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

import os, re, datetime
from . import dbreader, utils


def make_path_relative(filepath):
    if filepath[0] == '/':
        return filepath[1:]
    else:
        return filepath


def generate_soakdb_file(base_dir, input_dir):
    dbfile = os.path.join(base_dir, make_path_relative(input_dir), 'processing', 'database', 'soakDBDataFile.sqlite')
    return dbfile


def generate_xtal_dir(input_dir, xtal_name):
    """
    Generate the directory with the crystal data
    :param input_dir:
    :param xtal_name:
    :return:
    """
    xtal_dir = os.path.join(input_dir, 'processing', 'analysis', 'model_building', xtal_name)
    return xtal_dir


def prepend_base(base_dir, filepath):
    """
    Prepend the base path to the file path, if one exists.
    :param base_dir:
    :param filepath:
    :return:
    """
    if base_dir:
        full_inputpath = os.path.join(base_dir, make_path_relative(filepath))
    else:
        full_inputpath = filepath
    return full_inputpath


def generate_filenames(filepath, xtal_dir, output_dir):

    if os.path.isabs(filepath):
        # absolute file path
        inputpath = filepath
        outputpath = os.path.join(output_dir, make_path_relative(filepath))
    else:
        # relative path
        inputpath = os.path.join(xtal_dir, filepath)
        outputpath = os.path.join(output_dir, make_path_relative(xtal_dir), filepath)

    return inputpath, outputpath


class Validator:

    def __init__(self, base_dir, input_dirs, output_dir, target_name, logger=None):
        self.base_dir = base_dir
        self.input_dirs = input_dirs
        self.output_dir = output_dir
        self.target_name = target_name
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def validate_all(self):
        errors, warnings = self.validate_paths()
        if errors:
            self.logger.error('Path validation failed with', errors, 'errors')
            return None
        if warnings:
            self.logger.info('Continuing with', warnings, 'warnings')
        meta = self.validate_soakdb_data()
        return meta

    def validate_paths(self):
        """
        Validate that the inputs and outputs that are defined in the config file are valid.
        Warnings can probably be tolerated, but should be investigated.
        Errors should be fixed before trying again.

        :return: the number of errors and warnings
        """
        errors = 0
        warnings = 0
        if self.input_dirs:
            for input_dir in self.input_dirs:
                if self.base_dir:
                    input_path = self.base_dir + input_dir
                else:
                    input_path = input_dir

                if not os.path.exists(input_path):
                    self.logger.error('input_path does not exist:', input_path)
                    errors += 1
                elif not os.path.isdir(input_path):
                    self.logger.error('input_dir argument is not a directory:', input_path)
                    errors += 1
                else:
                    dbfile = generate_soakdb_file(self.base_dir, input_dir)
                    if not os.path.isfile(dbfile):
                        self.logger.error('SoakDB database not found:', dbfile)
                        errors += 1
        else:
            self.logger.error('input_dirs not defined in the config file')
            errors += 1

        if self.output_dir:
            if not os.path.exists(self.output_dir):
                self.logger.error('output_dir does not exist:', self.output_dir)
                errors += 1
            elif not os.path.isdir(self.output_dir):
                self.logger.error('output_dir argument is not a directory:', self.output_dir)
                errors += 1
        else:
            self.logger.error('Output dir not defined in the config file')
            errors += 1

        if not self.target_name:
            self.logger.error('target_name not defined')
            errors += 1
        else:
            if len(self.target_name) < 4:
                self.logger.error('target_name must have at least 4 characters:', self.target_name)
                errors += 1
            else:
                x = re.search("^[A-Za-z]+[A-Za-z0-9_\\-]*$", self.target_name)
                if not x:
                    self.logger.error('Invalid target_name:', self.target_name)
                    errors += 1

        self.logger.info('Path validation encountered {} errors and {} warnings'.format(errors, warnings))
        return errors, warnings

    def validate_soakdb_data(self):
        """
        Read info from the SoakDB database and verify that the necessary entries are present
        :return: The generated metadata
        """
        valid_ids = {}
        meta = {
            'run_on': str(datetime.datetime.now()),
            'input_dirs': self.input_dirs,
            'output_dir': self.output_dir,
            'crystals': valid_ids}

        for input_dir in self.input_dirs:
            dbfile = generate_soakdb_file(self.base_dir, input_dir)
            self.logger.info('Opening DB file:', dbfile)
            df = dbreader.filter_dbmeta(dbfile)
            count = 0
            processed = 0
            for index, row in df.iterrows():
                count += 1
                xtal_name = row['CrystalName']
                xtal_dir = generate_xtal_dir(input_dir, xtal_name)
                if not xtal_name:
                    self.logger.error('Crystal name not defined, cannot process row {}'.format(xtal_name))
                else:
                    self.logger.info('Processing crystal {} {}'.format(count, xtal_name))
                    missing_files = 0
                    expanded_files = []
                    for colname in ['RefinementPDB_latest', 'RefinementMTZ_latest', 'RefinementCIF']:
                        if not colname:
                            self.logger.error('File not defined for {}'.format(colname))
                            expanded_files.append(None)
                        else:
                            file = row[colname]
                            if file:
                                inputpath, outputpath = generate_filenames(file, xtal_dir, "")
                                full_inputpath = prepend_base(self.base_dir, inputpath)
                                ok = self._check_file_exists(full_inputpath)
                                if ok:
                                    expanded_files.append(inputpath)
                                else:
                                    expanded_files.append(None)
                                    missing_files += 1
                                    self.logger.warn('File {} for {} not found: {}'.format(colname, xtal_name, row[colname]))
                            else:
                                expanded_files.append(None)
                                self.logger.warn('Entry {} for {} not defined in SoakDB'.format(colname, xtal_name))

                    if missing_files > 0:
                        self.logger.warn('{} files for {} missing. Will not process'.format(missing_files, xtal_name))
                    else:
                        processed += 1
                        if xtal_name in valid_ids.keys():
                            self.logger.warn("Crystal {} already exists, it's data will be overriden".format(xtal_name))

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

    def _check_file_exists(self, filepath):
        """
        Check whether this path exists and is a file
        :param filepath:
        :return:
        """
        if not os.path.isfile(filepath):
            return False
        return True

