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

import sys
import argparse
import shutil
from pathlib import Path

import pandas as pd

from xchemalign import dbreader, collator, utils
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


class ManualCopier:
    def __init__(self, base_path: Path, input_path: Path, output_path: Path, excludes, logger):
        self.base_path = base_path
        self.input_path = input_path
        self.output_path = output_path
        self.excludes = excludes
        self.logger = logger
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
            self._log_error("base path must be defined")
        if not self.input_path:
            self._log_error("input path must be defined")
        elif self.input_path.is_absolute():
            self._log_error("input_path must be a path relative to base_path")
        return len(self.errors), len(self.warnings)

    def copy_files(self):
        if self.base_path and self.input_path.is_absolute():
            self.logger.warn("INFO: making input_path relative as a base_path is specified")
            self.input_path = self.input_path.relative_to("/")

        items = utils.collect_manual_files(self.base_path / self.input_path)
        self.logger.info("found {} manual inputs".format(len(items)))

        d = self.output_path / self.input_path
        d.mkdir(exist_ok=True, parents=True)

        for key, item in items.items():
            if key in self.excludes:
                self.logger.info("excluding manual crystal", key)
                continue
            pdb = item[0]
            mtz = item[1]
            cif = item[2]

            self.logger.info('copying manual files for', key)

            if pdb:
                shutil.copy2(pdb, d, follow_symlinks=True)
            if mtz:
                shutil.copy2(mtz, d, follow_symlinks=True)
            if cif:
                shutil.copy2(cif, d, follow_symlinks=True)


class FileCopier:
    def __init__(self, logger):
        self.logger = logger

    def do_copy(self, from_path, to_path):
        return shutil.copy2(from_path, to_path, follow_symlinks=True)

    def file_exists(self, path_to_check):
        return path_to_check.is_file()


class Copier:
    def __init__(
        self,
        base_path: Path,
        input_path: Path,
        output_path: Path,
        soakdb_file_path: Path,
        panddas_file_paths: list[Path],
        ref_datasets: list[str],
        logger=None,
    ):
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

        self.base_path = base_path
        self.input_path = input_path
        self.output_path = output_path
        self.soakdb_file_path = soakdb_file_path
        self.panddas_file_paths = panddas_file_paths
        self.ref_datasets = ref_datasets
        self.errors = []
        self.warnings = []
        self.copier = FileCopier(self.logger)

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def _log_warning(self, msg):
        self.logger.warn(msg)
        self.warnings.append(msg)

    def validate(self):
        if not self.base_path:
            self._log_error("base path must be defined")
        if not self.input_path:
            self._log_error("input path must be defined")
        elif self.input_path.is_absolute():
            self._log_error("input_path must be a path relative to base_path")
        if self.soakdb_file_path.is_absolute():
            self._log_error("soakdb_file_path must be a path relative to input_path")

        return len(self.errors), len(self.warnings)

    def check_path(self, path, expected_path):
        try:
            relp = path.relative_to(expected_path)
            return True
        except ValueError as ve:
            self.logger.warn('unexpected path for file:', path)
            return False

    def copy_files(self):
        if self.base_path and self.input_path.is_absolute():
            self.logger.warn("INFO: making input_path relative as a base_path is specified")
            self.input_path = self.input_path.relative_to("/")

        dbfile = self.base_path / self.input_path / self.soakdb_file_path
        dbfile_out = self.output_path / self.input_path / self.soakdb_file_path
        dbfile_out.parent.mkdir(exist_ok=True, parents=True)
        self.logger.info("copying soakdbfile", str(dbfile))
        f = self.copier.do_copy(dbfile, dbfile_out)
        if not f:
            self.logger.error("Failed to copy soakdb file {} to {}. Can't continue".format(dbfile, dbfile_out))
            sys.exit(1)

        self.logger.info("reading soakdb file", dbfile_out)
        df = dbreader.filter_dbmeta(dbfile_out, self.ref_datasets)
        self.logger.info("soakdb shape is", df.shape)
        count = 0
        num_files = 0
        num_csv = 0
        datasets = {}
        for index, row in df.iterrows():
            count += 1
            xtal_name = row["CrystalName"]
            status_str = str(row[Constants.SOAKDB_COL_REFINEMENT_OUTCOME])
            if status_str.startswith("7"):
                self.logger.info("ignoring {} as status is 7".format(xtal_name))
                continue

            xtal_dir_path = collator.generate_xtal_dir(self.input_path, xtal_name)
            self.logger.info("processing {} {}".format(count, xtal_name))
            expected_path = self.base_path / self.input_path / Constants.DEFAULT_MODEL_BUILDING_DIR

            file = row[Constants.SOAKDB_COL_PDB]
            result = self.copy_file_and_log(
                xtal_name, Constants.SOAKDB_COL_PDB, row[Constants.SOAKDB_COL_PDB], xtal_dir_path
            )
            num_files += result
            if result:
                datasets[xtal_name] = row[Constants.SOAKDB_COL_PDB]

            num_files += self.copy_file_and_log(
                xtal_name, Constants.SOAKDB_COL_MTZ_LATEST, row[Constants.SOAKDB_COL_MTZ_LATEST], xtal_dir_path
            )
            num_files += self.copy_file_and_log(
                xtal_name, Constants.SOAKDB_COL_MTZ_FREE, row[Constants.SOAKDB_COL_MTZ_FREE], xtal_dir_path
            )
            mmcif = row[Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST]
            if mmcif and mmcif != 'None':
                num_files += self.copy_file_and_log(
                    xtal_name,
                    Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST,
                    mmcif,
                    xtal_dir_path,
                )
            else:
                self._log_warning(
                    Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST + " not defined for crystal " + xtal_name
                )
            result = self.copy_file_and_log(
                xtal_name, Constants.SOAKDB_COL_CIF, row[Constants.SOAKDB_COL_CIF], xtal_dir_path
            )
            if result:
                # for the ligand CIF file also copy the corresponding PDB file
                filename = row[Constants.SOAKDB_COL_CIF]
                ok = self.copy_file(Path(filename).with_suffix(".pdb"), xtal_dir_path)
                if ok:
                    num_files += 1
                else:
                    self._log_warning("Ligand PDB file " + filename + " not copied for crystal " + xtal_name)

            # Copy files needed for PDB deposition if the status is 5 or 6
            # TODO - review which statuses should be copied
            if status_str.startswith("5") or status_str.startswith("6"):
                dp_log = row[Constants.SOAKDB_COL_DATA_PROCESSING_PATH_TO_LOGFILE]
                dp_prog = row[Constants.SOAKDB_COL_DATA_PROCESSING_PROGRAM]
                if dp_log and dp_log != 'None':
                    dp_log_p = Path(dp_log)
                    if dp_log_p.is_file():
                        self.logger.info('copying', str(dp_log))
                        num_files += self.copy_file_and_log(
                            xtal_name, Constants.SOAKDB_COL_DATA_PROCESSING_PATH_TO_LOGFILE, dp_log, xtal_dir_path
                        )

                    if dp_prog and dp_prog.lower() == 'autoproc' and dp_log_p.name.endswith('.log'):
                        # probably a broken symlink which should be to the aimless.log file
                        self.logger.info('copying aimless.log')
                        aimless = dp_log_p.parent / 'aimless.log'
                        if aimless.is_file():
                            num_files += self.copy_file_and_log(xtal_name, 'aimless.log', str(aimless), xtal_dir_path)

                    stats_cif_p = dp_log_p.parent / 'xia2.mmcif.bz2'
                    if stats_cif_p.is_file():
                        self.logger.info('copying', str(stats_cif_p))
                        num_files += self.copy_file_and_log(
                            xtal_name, 'xia2.mmcif.bz2', str(stats_cif_p), xtal_dir_path
                        )
                else:
                    self._log_warning("Data processing logfile not defined for crystal " + xtal_name)

        # copy the specified csv files with the panddas info
        self.logger.info("Copying", len(self.panddas_file_paths), "panddas csv files")
        copied_csv = []
        for panddas_path in self.panddas_file_paths:
            f = self.copy_csv(panddas_path)
            if f:
                copied_csv.append(f)
                num_csv += 1
            else:
                self.logger.error("Panddas CSV file not found", panddas_path)
                sys.exit(1)

        # copy the relevant panddas event map files
        num_ccp4 = self.copy_panddas(datasets, copied_csv)

        self.logger.info("Copied {} structure, {} csv, {} ccp4 files".format(num_files, num_csv, num_ccp4))

    def copy_file_and_log(self, xtal_name, col_name, col_value, xtal_dir_path):
        if col_value:
            path = Path(col_value)
            ok = self.copy_file(path, xtal_dir_path)
            if ok:
                return 1
            else:
                self._log_warning(
                    "File " + col_value + " not copied for column " + col_name + " and crystal " + xtal_name
                )
                return 0
        else:
            self._log_warning("Column value not defined for column " + col_name + " and crystal " + xtal_name)
            return 0

    def copy_file(self, filepath: Path, xtal_dir_path: Path):
        if filepath.is_absolute():
            inputpath_short = utils.make_path_relative(filepath)
            outputpath = self.output_path / filepath.relative_to("/")
        else:
            inputpath_short = xtal_dir_path / filepath
            outputpath = self.output_path / utils.make_path_relative(xtal_dir_path) / filepath

        inputpath_long = self.base_path / inputpath_short

        if not self.copier.file_exists(inputpath_long):
            self.logger.warn("file {} not found".format(inputpath_long))
            return False

        outputpath.parent.mkdir(exist_ok=True, parents=True)
        f = self.copier.do_copy(inputpath_long, outputpath)
        if not f:
            self.logger.warn("Failed to copy file {} to {}".format(inputpath_long, outputpath))
        return f

    def copy_csv(self, filepath: Path):
        """
        Copy a CSV files with the panddas data
        :param filepath:
        :return:
        """

        csv_src = self.base_path / self.input_path / filepath

        if not self.copier.file_exists(csv_src):
            self.logger.warn("File {} not found".format(csv_src))
            return False

        rel_path = utils.make_path_relative(self.input_path) / filepath
        dest_path = self.output_path / rel_path
        dest_dir_path = dest_path.parent
        dest_dir_path.mkdir(exist_ok=True, parents=True)
        dest_file = dest_dir_path / filepath.name
        self.logger.info("Copying CSV file {} to {}".format(csv_src, dest_file))
        f = self.copier.do_copy(csv_src, dest_file)
        if not f:
            self.logger.warn("Failed to copy CSV file {} to {}".format(csv_src, dest_file))
            return None
        return rel_path

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

        self.logger.info(len(datasets), "datasets")

        panddas_dict = {}
        ccp4_count = 0
        for panddas_file in panddas_files:
            self.logger.info("Handling", panddas_file)
            rel_pfp = self.input_path / panddas_file
            full_pfp = self.base_path / rel_pfp
            # if not full_pfp.exists():
            #     self.logger.warn("CSV file missing:", full_pfp)
            #     continue
            # self.logger.info("Reading CSV:", full_pfp)
            copied_panddas_path = self.output_path / panddas_file
            self.logger.info("Reading CSV", copied_panddas_path)
            df = pd.read_csv(copied_panddas_path)
            self.logger.info("Data frame shape:", df.shape)
            panddas_dict[panddas_file] = df

            for xtal_name, pdb_path in datasets.items():
                dataset_events = df[df[Constants.EVENT_TABLE_DTAG] == xtal_name]
                for idx, row in dataset_events.iterrows():
                    self.logger.info("event", idx, xtal_name)
                    event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
                    bdc = row[Constants.EVENT_TABLE_BDC]
                    found = False
                    event_map_paths = []
                    for template in Constants.EVENT_MAP_TEMPLATES:
                        event_map_path = (
                            panddas_file.parent.parent
                            / Constants.PROCESSED_DATASETS_DIR
                            / xtal_name
                            / template.format(dtag=xtal_name, event_idx=event_idx, bdc=bdc)
                        )
                        event_map_path_in = self.base_path / event_map_path
                        event_map_path_out = self.output_path / event_map_path
                        self.logger.info("EVENTMAP", panddas_file, event_map_path_in, event_map_path_out)

                        event_map_paths.append(event_map_path_out)
                        if self.copier.file_exists(event_map_path_in):
                            self.logger.info("PANDDAS", self.output_path, utils.make_path_relative(event_map_path_out))
                            event_map_path_out.parent.mkdir(exist_ok=True, parents=True)
                            f = self.copier.do_copy(event_map_path_in, event_map_path_out)
                            if not f:
                                self.logger.warn(
                                    "failed to copy event map file", event_map_path_in, event_map_path_out
                                )
                            else:
                                ccp4_count += 1
                                found = True
                    if not found:
                        self.logger.warn("event map file not found:", *event_map_paths)
        return ccp4_count

    def generate_file_paths(self, filepath: Path, xtal_dir_path: Path, output_path: Path):
        if filepath.is_absolute():
            inputpath = filepath
            outputpath = output_path / filepath.relative_to("/")
        else:
            inputpath = xtal_dir_path / filepath
            outputpath = output_path / utils.make_path_relative(xtal_dir_path) / filepath

        return inputpath, outputpath


def handle_inputs(base_dir, inputs, ref_datasets, output_dir, logger):
    """
    Works through the inputs and copies their data

    :param base_dir:
    :param inputs:
    :param ref_datasets:
    :param output_dir:
    :param logger:
    :return:
    """
    soakdb_files = []
    panddas_files = []
    input_dirs_model_building = []
    input_dirs_manual = []
    excludes_manual = []

    for input in inputs:
        logger.info("Looking at input", input.get(Constants.CONFIG_DIR))
        t = input.get(Constants.CONFIG_TYPE)
        if t == Constants.CONFIG_TYPE_MANUAL:
            input_dirs_manual.append(input.get(Constants.CONFIG_DIR))
            excludes_manual.append(input.get(Constants.CONFIG_EXCLUDE, []))
        else:
            input_dirs_model_building.append(input.get(Constants.CONFIG_DIR))
            soakdb_files.append(input.get(Constants.CONFIG_SOAKDB, 'processing/database/soakDBDataFile.sqlite'))
            panddas_files.append(input.get(utils.Constants.CONFIG_PANDDAS_EVENT_FILES, []))

    logger.info('Using model building input dirs:', input_dirs_model_building)
    logger.info('Using manual input dirs:', input_dirs_manual)

    copiers = []

    for i, input_dir in enumerate(input_dirs_model_building):
        msg = (
            "Running model_building copier. base_dir={}, input_dir={} output_dir={}, soakdbfile={}, panddas={}".format(
                base_dir, input_dir, output_dir, soakdb_files[i], ", ".join(panddas_files[i])
            )
        )
        logger.info(msg)

        c = Copier(
            Path(base_dir),
            Path(input_dir),
            Path(output_dir),
            Path(soakdb_files[i]),
            [Path(p) for p in panddas_files[i]],
            ref_datasets,
            logger=logger,
        )
        errors, warnings = c.validate()
        if errors:
            logger.error("There are errors, cannot continue")
            logger.report()
            logger.close()
            exit(1)
        else:
            copiers.append(c)

    for input_dir, excludes in zip(input_dirs_manual, excludes_manual):
        msg = "Running manual copier. base_dir={}, input_dir={} output_dir={}".format(base_dir, input_dir, output_dir)
        logger.info(msg)
        c = ManualCopier(Path(base_dir), Path(input_dir), Path(output_dir), excludes, logger)
        errors, warnings = c.validate()
        if errors:
            logger.error("There are errors, cannot continue")
            logger.report()
            logger.close()
            exit(1)
        else:
            copiers.append(c)

    for copier in copiers:
        logger.info("copying ...")
        copier.copy_files()


def main():
    # Example:
    #   python -m xchemalign.copier -c config.yaml -o inputs -l copier.log

    parser = argparse.ArgumentParser(description="copier")

    parser.add_argument("-c", "--config-file", required=True, default="config.yaml", help="Configuration file")
    parser.add_argument(
        "-b", "--base-dir", help="Base directory. If not specified then value from config.yaml is used"
    )
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("copier: ", args)
    utils.LOG = logger

    config = utils.read_config_file(args.config_file)
    ref_datasets = config.get(Constants.CONFIG_REF_DATASETS, [])

    if args.base_dir:
        base_dir = args.base_dir
    else:
        base_dir = config.get(Constants.CONFIG_BASE_DIR)

    inputs = config.get(Constants.CONFIG_INPUTS)
    if not inputs:
        logger.error("No inputs defined in config file")
        sys.exit(1)

    handle_inputs(base_dir, inputs, ref_datasets, args.output_dir, logger)

    logger.report()
    logger.close()


if __name__ == "__main__":
    main()
