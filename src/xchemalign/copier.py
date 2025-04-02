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

from paramiko import SSHClient
from scp import SCPClient

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


class FileCopier:
    def __init__(self, logger):
        self.logger = logger

    def do_copy(self, from_path, to_path):
        return shutil.copy2(from_path, to_path, follow_symlinks=True)

    def file_exists(self, path_to_check):
        return path_to_check.is_file()


class SCPCopier:
    def __init__(self, logger, username, server='ssh.diamond.ac.uk', key_filename=None):
        self.logger = logger
        self.server = server
        self.ssh = SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.connect(self.server, username=username, key_filename=key_filename)
        self.scp = SCPClient(self.ssh.get_transport())

    def __del__(self):
        try:
            self.scp.close()
        except NameError:
            pass

    def do_copy(self, from_path, to_path):
        try:
            self.logger.info("Copying", from_path, "to", to_path)
            self.scp.get(str(from_path), str(to_path))
            return to_path.is_file()
        except IOError:
            # file not found?
            return False

    def file_exists(self, path_to_check):
        try:
            sftp = self.ssh.open_sftp()
            sftp.chdir(str(path_to_check.parent))  # changing to directory we want to search
            sftp.stat(str(path_to_check))  # checking whether file exists or not
            return True
        except IOError:
            return False


class Copier:
    def __init__(
        self,
        base_path: Path,
        input_path: Path,
        output_path: Path,
        soakdb_file_path: Path,
        panddas_file_paths: list[Path],
        mode: str,
        ref_datasets: list[str],
        scp_user: str = None,
        scp_server="ssh.diamond.ac.uk",
        scp_key=None,
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
        self.mode = mode
        self.ref_datasets = ref_datasets
        self.errors = []
        self.warnings = []
        if mode == 'copy':
            self.copier = FileCopier(self.logger)
        elif mode == 'scp':
            self.copier = SCPCopier(self.logger, scp_user, server=scp_server, key_filename=scp_key)
        else:
            raise ValueError("Invalid copy mode:", mode)

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
        f = self.copier.do_copy(dbfile, dbfile_out)
        if not f:
            self.logger.error("Failed to copy soakdb file {} to {}. Can't continue".format(dbfile, dbfile_out))
            sys.exit(1)

        self.logger.info("reading soakdb file", dbfile_out)
        df = dbreader.filter_dbmeta(dbfile_out, self.ref_datasets)
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
            if file:
                path = Path(file)
                self.check_path(path, expected_path)
                ok = self.copy_file(path, xtal_dir_path)
                if ok:
                    num_files += 1
                    datasets[xtal_name] = path
                    # if PDB is OK then continue with the other files
                    file = row[Constants.SOAKDB_COL_MTZ]
                    if file:
                        path = Path(file)
                        ok = self.copy_file(path, xtal_dir_path)
                        if ok:
                            num_files += 1
                    file = row[Constants.SOAKDB_COL_CIF]
                    if file:
                        path = Path(file)
                        ok = self.copy_file(path, xtal_dir_path)
                        if ok:
                            num_files += 1
                            # for the ligand CIF file also copy the corresponding PDB file
                            ok = self.copy_file(path.with_suffix(".pdb"), xtal_dir_path)
                            if ok:
                                num_files += 1

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


def main():
    parser = argparse.ArgumentParser(description="copier")

    parser.add_argument("-c", "--config-file", default="config.yaml", help="Configuration file")

    parser.add_argument("-b", "--base-dir", help="Base directory. If running against the Diamond file system use /")
    parser.add_argument(
        "-i",
        "--input-dir",
        help="Input directory (relative to base-dir) e.g. the dir with the data for your visit. "
        + "e.g. dls/labxchem/data/2020/lb18145-153",
    )
    parser.add_argument(
        "-s",
        "--soakdb-file",
        help="Path to soakdb file relative to input-dir. Default is processing/database/soakDBDataFile.sqlite",
    )
    parser.add_argument(
        "-p", "--panddas-files", nargs="*", help="Paths to CSV files with panddas data relative to input-dir"
    )
    parser.add_argument("-m", "--mode", required=True, choices=['copy', 'scp'], help="Mode of file copying")
    parser.add_argument("--scp-username", help="SCP username")
    parser.add_argument("--scp-server", default="ssh.diamond.ac.uk", help="SCP server")
    parser.add_argument("--scp-key", help="SSH key")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-r", "--ref-datasets", nargs="*", help="Names of any reference datasets")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("copier: ", args)
    utils.LOG = logger

    input_dir = args.input_dir
    base_dir = None
    soakdbfiles = []
    panddas_files = []
    input_dirs = []
    scp_server = None
    scp_username = None
    scp_key = None
    ref_datasets = []
    if args.config_file:
        config = utils.read_config_file(args.config_file)
        scp_config = config.get('scp')
        if scp_config:
            scp_server = scp_config.get('server', 'ssh.diamond.ac.uk')
            scp_username = scp_config.get('username')
            scp_key = scp_config.get('key')
            base_dir = scp_config.get('base_dir', '/')

        ref_datasets = config.get(Constants.CONFIG_REF_DATASETS, [])

        inputs = config.get('inputs')
        if not inputs:
            logger.error("No inputs defined in config file")
            sys.exit(1)

        for input in inputs:
            logger.info("Looking at input", input.get('dir'))
            if not input_dir or input.get('dir') == input_dir:
                t = input.get('type')
                if t != 'model_building':
                    logger.warn("Only copying of model_building types is currently supported. You specified type", t)
                    continue
                logger.info("Adding input", input.get('dir'))
                input_dirs.append(input.get('dir'))
                soakdbfiles.append(input.get('soakdb', 'processing/database/soakDBDataFile.sqlite'))
                panddas_files.append(input.get(utils.Constants.CONFIG_PANDDAS_EVENT_FILES, []))

        # check we have at least one input
        if len(input_dirs) == 0:
            if input_dir:
                # a specific input was requested but was not found
                logger.error("Input {} not defined in config file".format(input_dir))
            else:
                # no input was requested and none were found
                logger.error("No inputs found in config file")
            sys.exit(1)

    if args.base_dir:
        base_dir = args.base_dir
    if base_dir is None:
        logger.error("base_dir must either be specified in the scp section of the config file or as an argument")
        sys.exit(1)

    # CLI soakdb_file arg override what is in the config file but only if a single input is defined
    if args.soakdb_file:
        if len(input_dirs) != 1:
            logger.error(
                "soakdb_file command line argument is ambiguous when processing multiple inputs. "
                + "Specify this in each input section of the config file instead."
            )
            sys.exit(1)
        soakdbfiles[0] = args.soakdb_file

    # CLI panddas_files arg override what is in the config file but only if a single input is defined
    if args.panddas_files:
        if len(input_dirs) != 1:
            logger.error(
                "panddas_files command line argument is ambiguous when processing multiple inputs. "
                + "Specify this in each input section of the config file instead."
            )
            sys.exit(1)
        panddas_files[0] = args.panddas_files

    if args.ref_datasets:
        ref_datasets = args.ref_datasets

    # CLI scp args override what is in the config file
    if args.scp_server:
        scp_server = args.scp_server
    if args.scp_username:
        scp_username = args.scp_username
    if args.scp_key:
        scp_key = args.scp_key

    for i, input_dir in enumerate(input_dirs):
        msg = (
            "Running copier using mode {}. base_dir={}, input_dir={} output_dir={}, soakdbfile={}, panddas={}".format(
                args.mode, base_dir, input_dir, args.output_dir, soakdbfiles[i], ", ".join(panddas_files[i])
            )
        )
        logger.info(msg)

        c = Copier(
            Path(base_dir),
            Path(input_dir),
            Path(args.output_dir),
            Path(soakdbfiles[i]),
            [Path(p) for p in panddas_files[i]],
            args.mode,
            ref_datasets,
            scp_server=scp_server,
            scp_user=scp_username,
            scp_key=scp_key,
            logger=logger,
        )
        errors, warnings = c.validate()
        if errors:
            logger.error("There are errors, cannot continue")
            logger.report()
            logger.close()
            ex = 1
        else:
            c.copy_files()
            ex = 0

        logger.report()
        logger.close()
        exit(ex)


if __name__ == "__main__":
    main()
