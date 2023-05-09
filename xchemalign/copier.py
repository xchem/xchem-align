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

import os, argparse, shutil
from pathlib import Path

import pandas as pd

from . import dbreader, processor, utils
from .utils import Constants



def _generate_path(base_path: Path, input_path: Path, file_path):
    if input_path:
        if base_path:
            return base_path / input_path / file_path
        else:
            return input_path / file_path
    else:
        if base_path:
            return base_path / file_path
        else:
            return file_path

class Copier(processor.Processor):

    def __init__(self, base_path: Path, input_path:Path, output_path: Path, soakdb_file_path: Path, panddas_file_paths: list[Path], logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

        self.base_path = base_path
        self.input_path = input_path
        self.output_path = output_path
        self.soakdb_file_path = soakdb_file_path
        self.panddas_file_paths = panddas_file_paths
        self.errors = []
        self.warnings = []

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def _log_warning(self, msg):
        self.logger.warn(msg)
        self.warnings.append(msg)
    def validate(self):
        if not self.base_path:
            self._log_error('base path must be defined')
        if not self.input_path:
            self._log_error('input path must be defined')
        elif self.input_path.is_absolute():
            self._log_error('input_path must be a path relative to base_path')
        if self.soakdb_file_path.is_absolute():
            self._log_error('soakdb_file_path must be a path relative to input_path')

        return len(self.errors), len(self.warnings)

    def copy_files(self):

        if self.base_path and self.input_path.is_absolute():
            self.logger.warn('INFO: making input_path relative as a base_path is specified')
            self.input_path = self.input_path.relative_to('/')

        dbfile = self.base_path / self.input_path / self.soakdb_file_path
        dbfile_out = self.output_path / self.input_path / self.soakdb_file_path
        dbfile_out.parent.mkdir(exist_ok=True, parents=True)
        f = shutil.copy2(dbfile, dbfile_out, follow_symlinks=True)

        self.logger.info('reading soakdb file', dbfile)
        df = dbreader.filter_dbmeta(dbfile)
        count = 0
        num_files = 0
        num_csv = 0
        datasets = {}
        for index, row in df.iterrows():
            count += 1
            xtal_name = row['CrystalName']
            xtal_dir_path = processor.generate_xtal_dir(self.input_path, xtal_name)
            self.logger.info('processing {} {}'.format(count, xtal_name))

            file = row['RefinementPDB_latest']
            if file:
                path = Path(file)
                ok = self.copy_file(path, xtal_dir_path)
                if ok:
                    num_files += 1
                    datasets[xtal_name] = path
                    # if PDB is OK then continue with the other files
                    file = row['RefinementMTZ_latest']
                    if file:
                        path = Path(file)
                        ok = self.copy_file(path, xtal_dir_path)
                        if ok:
                            num_files += 1
                    file = row['RefinementCIF']
                    if file:
                        path = Path(file)
                        ok = self.copy_file(path, xtal_dir_path)
                        if ok:
                            num_files += 1
                            # for the ligand CIF file also copy the corresponding PDB file
                            ok = self.copy_file(path.with_suffix('.pdb'), xtal_dir_path)
                            if ok:
                                num_files += 1

        # copy the specified csv files with the panddas info
        self.logger.info('Copying panddas csv files')
        for panddas_path in self.panddas_file_paths:
            ok = self.copy_csv(panddas_path)
            if ok:
                num_csv += 1

        # copy the relevant panddas event map files
        num_ccp4 = self.copy_panddas(datasets, self.panddas_file_paths)

        self.logger.info('Copied {} structure, {} csv, {} ccp4 files'.format(num_files, num_csv, num_ccp4))

    def copy_file(self, filepath: Path, xtal_dir_path: Path):

        if filepath.is_absolute():
            inputpath_short =  utils.make_path_relative(filepath)
            outputpath = self.output_path / filepath.relative_to('/')
        else:
            inputpath_short = xtal_dir_path / filepath
            outputpath = self.output_path / utils.make_path_relative(xtal_dir_path) / filepath

        inputpath_long = self.base_path / inputpath_short

        # print('copying', inputpath_long, outputpath)

        if not inputpath_long.is_file():
            self.logger.warn('file {} not found'.format(inputpath_long))
            return False

        outputpath.parent.mkdir(exist_ok=True, parents=True)
        f = shutil.copy2(inputpath_long, outputpath, follow_symlinks=True)
        if not f:
            self.logger.warn('Failed to copy file {} to {}'.format(inputpath_long, outputpath))
            return False

        return True

    def copy_csv(self, filepath: Path):
        """
        Copy a CSV files with the panddas data
        :param filepath:
        :return:
        """

        csv_src = self.base_path / self.input_path / filepath

        if not csv_src.is_file():
            self.logger.warn('File {} not found'.format(csv_src))
            return False

        dest_path = self.output_path / utils.make_path_relative(self.input_path) / filepath
        dest_dir_path = dest_path.parent
        dest_dir_path.mkdir(exist_ok=True, parents=True)
        f = shutil.copy2(csv_src, dest_dir_path, follow_symlinks=True)
        if not f:
            self.logger.warn('Failed to copy file {} to {}'.format(csv_src, dest_dir_path))
            return False
        return True

    def copy_panddas(self, datasets: dict[str, Path], panddas_files: list[Path]):
        """
        Find the PanDDAs event maps (.ccp4 files) that are relevant to the data contained in the CSV files
        can copy then to the corresponding locations in the output dir.
        :param datasets: Dict, keys being the crystal name (e.g. Mpro-x1234) and values being the path to the PDB file
        :param panddas_files: A list of paths to the CSV files with the panddas data
        :return: The number of .ccp4 files copies
        """

        if not panddas_files:
            return

        self.logger.info(len(datasets), 'datasets')

        panddas_dict = {}
        ccp4_count = 0
        for panddas_file in panddas_files:
            pfp = self.base_path / self.input_path / panddas_file
            self.logger.info('Reading CSV:', pfp)
            df = pd.read_csv(pfp)
            self.logger.info('Data frame shape:', df.shape)
            panddas_dict[pfp] = df

            for xtal_name, pdb_path in datasets.items():

                dataset_events = df[df[Constants.EVENT_TABLE_DTAG] == xtal_name]
                for idx, row in dataset_events.iterrows():
                    self.logger.info('event', idx, xtal_name)
                    event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
                    bdc = row[Constants.EVENT_TABLE_BDC]
                    found = False
                    event_map_paths = []
                    for template in Constants.EVENT_MAP_TEMPLATES:
                        event_map_path = pfp.parent.parent / Constants.PROCESSED_DATASETS_DIR / xtal_name / template.format(
                            dtag=xtal_name,
                            event_idx=event_idx,
                            bdc=bdc
                        )
                        event_map_paths.append(event_map_path)
                        if event_map_path.exists() and event_map_path.is_file():
                            outfile = self.output_path / utils.make_path_relative(self.input_path) / event_map_path
                            outfile.parent.mkdir(exist_ok=True, parents=True)
                            f = shutil.copy2(event_map_path, outfile, follow_symlinks=True)
                            if not f:
                                self.logger.warn('failed to copy event map file', event_map_path)
                            else:
                                ccp4_count += 1
                                found = True
                    if not found:
                        self.logger.warn('event map file not found:', *event_map_paths)
        return ccp4_count

    def generate_file_paths(self, filepath: Path, xtal_dir_path: Path, output_path: Path):

        if filepath.is_absolute():
            inputpath = filepath
            outputpath = output_path / filepath.relative_to('/')
        else:
            inputpath = xtal_dir_path / filepath
            outputpath = output_path / xtal_dir_path / filepath

        return inputpath, outputpath


def main():
    parser = argparse.ArgumentParser(description='copier')

    parser.add_argument('-b', '--base-dir', required=True, help="Base directory")
    parser.add_argument('-i', '--input-dir', required=True, help="Input directory (relative to base-dir)")
    parser.add_argument('-s', '--soakdb-file', default='processing/database/soakDBDataFile.sqlite',
                        help="Path to soakdb file relative to input-dir")
    parser.add_argument('-p', '--panddas-files', nargs='*', help="Path to CSV files with panddas data relative to input-dir")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("copier: ", args)

    if args.panddas_files:
        panddas_paths = [Path(p) for p in args.panddas_files]
    else:
        panddas_paths = []

    if args.base_dir:
        base_path = Path(args.base_dir)
    else:
        base_path = None

    c = Copier(base_path, Path(args.input_dir), Path(args.output_dir), Path(args.soakdb_file), panddas_paths)
    errors, warnings = c.validate()
    if errors:
        print('There are errors, cannot continue')
        exit(1)
    else:
        c.copy_files()


if __name__ == "__main__":
    main()
