import os, re, datetime, collections
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
                print('input_dir does not exist:', self.input_dir)
                errors += 1
            elif not os.path.isdir(self.input_dir):
                print('input_dir argument is not a directory:', self.input_dir)
                errors += 1
        else:
            print('Input dir not defined. Use the --input-dir to define this')
            errors += 1

        if self.output_dir:
            if not os.path.exists(self.output_dir):
                print('output_dir does not exist:', self.output_dir)
                errors += 1
            elif not os.path.isdir(self.output_dir):
                print('output_dir argument is not a directory:', self.output_dir)
                errors += 1
        else:
            print('Output dir not defined. Use the --output-dir to define this')
            errors += 1

        if not self.target_name:
            print('target_name not defined')
            errors += 1
        else:
            if len(self.target_name) < 4:
                print('target_name must have at least 4 characters:', self.target_name)
                errors += 1
            else:
                x = re.search("^[A-Za-z]+[A-Za-z0-9_\\-]*$", self.target_name)
                if not x:
                    print('Invalid target_name:', self.target_name)
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
        for index, row in df.iterrows():
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
                            'ligand_cif': row['RefinementCIF'],}}
        return meta

    def _check_file_exists(self, row, colname):
        filepath = row[colname]
        if filepath:
            if self.input_dir:
                path = self.input_dir + filepath
            else:
                path = filepath
            if not os.path.isfile(path):
                #self.logger.log('{} file {} not found'.format(colname, filepath), level=1)
                return False
        else:
            self.logger.log('{} file not defined'.format(colname), level=1)
            return False
        return True


def main():

    v = Validator('data', '', '', 'data/dls/labxchem/data/lb18145/lb18145-216/processing/database/soakDBDataFile.sqlite')

    meta = v.validate_metadata()
    print(meta)


if __name__ == "__main__":
    main()