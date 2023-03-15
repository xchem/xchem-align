import os, re, datetime, argparse, shutil
from . import dbreader
from . import utils


def generate_soakdb_file(base_dir, input_dir):
    dbfile = base_dir + '/' + input_dir + '/processing/database/soakDBDataFile.sqlite'
    return dbfile


def generate_xtal_dir(input_dir, xtal_name):
    xtal_dir = os.path.join(input_dir, 'processing', 'analysis', 'model_building', xtal_name)
    return xtal_dir


def prepend_base(base_dir, filepath):
    if base_dir:
        full_inputpath = base_dir + '/' + filepath
    else:
        full_inputpath = filepath
    return full_inputpath


def generate_filenames(filepath, xtal_dir, output_dir):

    if filepath[0] == '/':
        # absolute file path
        inputpath = filepath
        outputpath = output_dir + filepath
    else:
        # relative path
        inputpath = xtal_dir + '/' + filepath
        outputpath = output_dir + '/' + xtal_dir + '/' + filepath

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
            self.logger.log('Path validation failed with', errors, 'errors', level=2)
            return None
        if warnings:
            self.logger.log('Continuing with', warnings, 'warnings', level=0)
        meta = self.validate_metadata()
        return meta

    def validate_paths(self):
        errors = 0
        warnings = 0
        if self.input_dirs:
            for input_dir in self.input_dirs:
                if self.base_dir:
                    input_path = self.base_dir + input_dir
                else:
                    input_path = input_dir

                if not os.path.exists(input_path):
                    self.logger.log('input_path does not exist:', input_path, level=2)
                    errors += 1
                elif not os.path.isdir(input_path):
                    self.logger.log('input_dir argument is not a directory:', input_path, level=2)
                    errors += 1
                else:
                    dbfile = generate_soakdb_file(self.base_dir, input_dir)
                    if not os.path.isfile(dbfile):
                        self.logger.log('SoakDB database not found:', dbfile, level=2)
                        errors += 1

        else:
            self.logger.log('input_dirs not defined in the config file', level=2)
            errors += 1

        if self.output_dir:
            if not os.path.exists(self.output_dir):
                self.logger.log('output_dir does not exist:', self.output_dir, level=2)
                errors += 1
            elif not os.path.isdir(self.output_dir):
                self.logger.log('output_dir argument is not a directory:', self.output_dir, level=2)
                errors += 1
        else:
            self.logger.log('Output dir not defined in the config file', level=2)
            errors += 1

        if not self.target_name:
            self.logger.log('target_name not defined', level=2)
            errors += 1
        else:
            if len(self.target_name) < 4:
                self.logger.log('target_name must have at least 4 characters:', self.target_name, level=2)
                errors += 1
            else:
                x = re.search("^[A-Za-z]+[A-Za-z0-9_\\-]*$", self.target_name)
                if not x:
                    self.logger.log('Invalid target_name:', self.target_name, level=2)
                    errors += 1

        self.logger.log('Path validation encounters {} errors and {} warnings'.format(errors, warnings), level=0)
        return errors, warnings

    def validate_metadata(self):
        meta = {'run_on': str(datetime.datetime.now()), 'input_dirs': self.input_dirs, 'output_dir': self.output_dir}
        for input_dir in self.input_dirs:
            dbfile = generate_soakdb_file(self.base_dir, input_dir)
            print('Opening DB file:', dbfile)
            df = dbreader.filter_dbmeta(dbfile)
            valid_ids = {}
            meta['crystals'] = valid_ids
            count = 0
            processed = 0
            for index, row in df.iterrows():
                count += 1
                xtal_name = row['CrystalName']
                xtal_dir = generate_xtal_dir(input_dir, xtal_name)
                if not xtal_name:
                    self.logger.log('Crystal name not defined, cannot process row {}'.format(count), level=2)
                else:
                    missing_files = 0
                    expanded_files = []
                    for colname in ['RefinementPDB_latest', 'RefinementMTZ_latest', 'RefinementCIF']:
                        if not colname:
                            self.logger.log('File not defined for {}'.format(colname), level=2)
                            expanded_files.append(None)
                        else:
                            file = row[colname]
                            inputpath, outputpath = generate_filenames(file, xtal_dir, "")
                            full_inputpath = prepend_base(self.base_dir, inputpath)
                            ok = self._check_file_exists(full_inputpath)
                            if ok:
                                expanded_files.append(inputpath)
                            else:
                                expanded_files.append(None)
                                missing_files += 1
                                self.logger.log('File for {} not found: {}'.format(colname, row[colname]), level=1)

                    if missing_files > 0:
                        self.logger.log('{} files for {} missing. Will not process'.format(missing_files, xtal_name),
                                        level=1)
                    else:
                        processed += 1
                        if xtal_name in valid_ids.keys():
                            self.logger.log("Crystal {} already exists, it's data will be overriden".format(xtal_name),
                                            level=1)

                        self.logger.log('Adding crystal', xtal_name, level=0)
                        data = {}
                        valid_ids[xtal_name] = data
                        last_updated = row['LastUpdated']
                        if last_updated:
                            data['last_updated'] = last_updated
                        data['crystallographic_files'] = {
                            'xtal_pdb': expanded_files[0],
                            'xtal_mtz': expanded_files[1],
                            'ligand_cif': expanded_files[2]}

            self.logger.log('Handled {} rows from database, {} were valid'.format(count, processed), level=0)

        return meta

    def _check_file_exists(self, filepath):
        if not os.path.isfile(filepath):
            self.logger.log('File {} not found'.format(filepath), level=1)
            return False
        return True


def main():

    parser = argparse.ArgumentParser(description='processor')

    parser.add_argument('-i', '--input-dir', required=True, help="Input directory")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")
    parser.add_argument('-d', '--db-file', required=True, help="Sqlite DB file")

    args = parser.parse_args()
    print("validator: ", args)

    # v = Validator([args.input_dir], args.output_dir, 'targetname')
    #
    # v.copy_files(args.db_file)


if __name__ == "__main__":
    main()
