import os, re, datetime, argparse, shutil
from . import dbreader
from . import utils

class Validator:

    def __init__(self, input_dir, output_dir, target_name, dbfile, logger=None):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.target_name = target_name
        self.dbfile = dbfile
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def copy_files(self):
        print('copy files')
        df = dbreader.filter_dbmeta(self.dbfile)
        count = 0
        for index, row in df.iterrows():
            count += 1
            print('processing {} {}'.format(count, row['CrystalName']))
            for colname in ['RefinementCIF', 'RefinementPDB_latest', 'RefinementMTZ_latest']:
                ok = self._copy_file(row, colname)
                if ok:
                    print('copy ok')

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
        if self.input_dir:
            if not os.path.exists(self.input_dir):
                self.logger.log('input_dir does not exist:', self.input_dir, level=2)
                errors += 1
            elif not os.path.isdir(self.input_dir):
                self.logger.log('input_dir argument is not a directory:', self.input_dir, level=2)
                errors += 1
        else:
            self.logger.log('Input dir not defined. Use the --input-dir to define this', level=2)
            errors += 1

        if self.output_dir:
            if not os.path.exists(self.output_dir):
                self.logger.log('output_dir does not exist:', self.output_dir, level=2)
                errors += 1
            elif not os.path.isdir(self.output_dir):
                self.logger.log('output_dir argument is not a directory:', self.output_dir, level=2)
                errors += 1
        else:
            self.logger.log('Output dir not defined. Use the --output-dir to define this', level=2)
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
        meta = {}
        meta['run_on'] = str(datetime.datetime.now())
        meta['input_dir'] = self.input_dir
        meta['output_dir'] = self.output_dir
        df = dbreader.filter_dbmeta(self.dbfile)
        valid_ids = {}
        meta['crystals'] = valid_ids
        count = 0
        for index, row in df.iterrows():
            count += 1
            for colname in ['RefinementCIF', 'RefinementPDB_latest', 'RefinementMTZ_latest']:
                ok = self._check_file_exists(row, colname)
                if ok:
                    xtal_name = row['CrystalName']
                    if not xtal_name:
                        self.logger.log('Crystal name not found', level=1)
                    else:
                        valid_ids[xtal_name] = {'crystallographic_files': {
                            'xtal_pdb': row['RefinementPDB_latest'],
                            'xtal_mtz': row['RefinementMTZ_latest'],
                            'ligand_cif': row['RefinementCIF']}}
        self.logger.log('Handled {} rows from database'.format(count), level=1)
        return meta

    def _check_file_exists(self, row, colname):
        filepath = row[colname]
        if filepath:
            if self.input_dir:
                path = self.input_dir + filepath
            else:
                path = filepath
            if not os.path.isfile(path):
                self.logger.log('{} file {} not found'.format(colname, filepath), level=1)
                return False
        else:
            self.logger.log('{} file not defined'.format(colname), level=1)
            return False
        return True

    def _copy_file(self, row, colname):
        filepath = row[colname]
        if filepath:
            if self.input_dir:
                path = self.input_dir + filepath
            else:
                path = filepath
            if not os.path.isfile(path):
                self.logger.log('{} file {} not found'.format(colname, path), level=1)
                return False
        else:
            self.logger.log('{} file not defined'.format(colname), level=1)
            return False
        outputpath = self.output_dir + filepath
        os.makedirs(os.path.dirname(outputpath), exist_ok=True)
        f = shutil.copy2(path, outputpath, follow_symlinks=True)
        if not f:
            self.logger.log('Failed to copy CIF file {} to {}'.format(path, outputpath), level=2)
            return False
        return True


def main():

    parser = argparse.ArgumentParser(description='processor')

    parser.add_argument('-i', '--input-dir', required=True, help="Input directory")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")
    parser.add_argument('-d', '--db-file', required=True, help="Sqlite DB file")

    args = parser.parse_args()
    print("validator: ", args)

    v = Validator(args.input_dir, args.output_dir, '', args.db_file)

    v.copy_files()


if __name__ == "__main__":
    main()
