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

import sys, os, atexit, hashlib
from pathlib import Path
import datetime
import yaml, json

_DATETIME_FORMAT = '%Y-%m-%d %H:%M:%S'

class Constants:
    EVENT_TABLE_DTAG = "dtag"
    EVENT_TABLE_EVENT_IDX = "event_idx"
    EVENT_TABLE_X = "x"
    EVENT_TABLE_Y = "y"
    EVENT_TABLE_Z = "z"
    EVENT_TABLE_BDC = "1-BDC"
    LIGAND_NAMES = ["LIG", "XXX"]
    PROCESSED_DATASETS_DIR = "processed_datasets"
    EVENT_MAP_TEMPLATES = [
        "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.ccp4",
        "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.native.ccp4"
    ]
    XTALFORMS_FILENAME = 'xtalforms.yaml'
    METADATA_XTAL_FILENAME = 'meta_collator.yaml'
    METADATA_ALIGN_FILENAME = 'meta_aligner.yaml'
    VERSION_DIR_PREFIX = 'upload_'
    DEFAULT_SOAKDB_PATH = 'processing/database/soakDBDataFile.sqlite'
    DEFAULT_MODEL_BUILDING_DIR = 'processing/analysis/model_building'
    CONFIG_INPUTS = 'inputs'
    CONFIG_TYPE = 'type'
    CONFIG_DIR = 'dir'
    CONFIG_SOAKDB = 'soakdb'
    CONFIG_TYPE_MODEL_BUILDING = 'model_building'
    CONFIG_TYPE_MANUAL = 'manual'
    CONFIG_BASE_DIR = 'base_dir'
    CONFIG_OUTPUT_DIR = 'output_dir'
    CONFIG_TARGET_NAME = 'target_name'
    META_RUN_ON = 'run_on'
    META_INPUT_DIRS = 'input_dirs'
    META_LAST_UPDATED = 'last_updated'
    META_STATUS = 'status'
    META_STATUS_SUPERSEDES = 'supersedes'
    META_STATUS_UNCHANGED = 'unchanged'
    META_STATUS_NEW = 'new'
    META_STATUS_DEPRECATED = 'deprecated'
    META_REASON = 'reason'
    META_XTALS = 'crystals'
    META_FILE = 'file'
    META_SHA256 = 'sha256'
    META_XTAL_FILES = 'crystallographic_files'
    META_ALIGNED_FILES = 'aligned_files'
    META_XTAL_PDB = 'xtal_pdb'
    META_XTAL_MTZ = 'xtal_mtz'
    META_XTAL_CIF = 'ligand_cif'
    META_BINDING_EVENT = 'panddas_event_files'
    META_PROT_MODEL = 'model'
    META_PROT_CHAIN = 'chain'
    META_PROT_RES = 'res'
    META_PROT_INDEX = 'index'
    META_PROT_BDC = 'bdc'
    META_AIGNED_STRUCTURE = 'structure'
    META_AIGNED_ARTEFACTS = 'artefacts'
    META_AIGNED_EVENT_MAP = 'event_map'
    META_AIGNED_X_MAP = 'x_map'
    META_CONFORMER_SITES = 'conformer_sites'
    META_CONFORMER_SITE_NAME = 'name'
    META_CONFORMER_SITE_REFERENCE_LIG = 'lig_ref'
    META_CONFORMER_SITE_RESIDUES = 'residues'
    META_CONFORMER_SITE_MEMBERS = 'members'
    META_CANONICAL_SITES = 'canon_sites'
    META_CANONICAL_SITE_REF_SUBSITE = 'canon_site_ref_site'
    META_CANONICAL_SITE_CONFORMER_SITES = 'canon_site_conf_sites'
    META_CANONICAL_SITE_RESIDUES = 'site_residues'
    META_CANONICAL_SITE_MEMBERS = 'site_members'
    META_XTALFORM_SITES = 'xtalform_sites'
    META_XTALFORM_SITE_XTALFORM_ID = 'xtalform_id'
    META_XTALFORM_SITE_CANONICAL_SITE_ID = 'xtalform_site_id'
    META_XTALFORM_SITE_LIGAND_CHAIN = 'lig_chain'
    META_XTALFORM_SITE_MEMBERS = 'xtalform_members'
    META_ASSIGNED_XTALFORM = 'assigned_xtalform'
    META_CHAIN = "chain"
    META_RESIDUE = "residue"
    META_DTAG = "dtag"
    META_XTALFORM_REFERENCE = "xtalform_ref"
    META_XTALFORM_SPACEGROUP = "xtalform_space_group"
    META_XTALFORM_CELL = "xtalform_cell"
    META_ASSEMBLIES_XTALFORMS = "assemblies_xtalforms"
    META_XTALFORMS = "xtalforms"
    META_ASSEMBLIES = 'assemblies'
    SOAKDB_XTAL_NAME = 'CrystalName'
    SOAKDB_COL_PDB = 'RefinementPDB_latest'
    SOAKDB_COL_MTZ = 'RefinementMTZ_latest'
    SOAKDB_COL_CIF = 'RefinementCIF'
    SOAKDB_COL_LAST_UPDATED = 'LastUpdatedDate'

class Logger:
    """
    Logger class that allows to write lines to the console and/or a file.
    """

    def __init__(self, logfile=None, console=sys.stderr, level=0):
        """

        :param logfilename: The name of a file to log to. If none then messages are not written to a file
        :param console: Whether to write messages to the console. The default is to write to sys.stderr, but you can
        specify sys.stdout or None instead.
        :param level: What types of message to log. 0 = everything, 1 = WARNING and ERROR, 2 = ERROR only
        """
        self.console = console
        self.level = 0
        self.infos = 0
        self.warnings = 0
        self.errors = 0
        if logfile:
            self.logfile = open(logfile, 'w')
            self.closed = False
        else:
            self.logfile = None
            self.closed = True
        atexit.register(self.close)
        x = datetime.datetime.now()
        self.log('Initialising logging at level {} at {}'.format(level, x), level=0)
        self.level = level

    def close(self):
        if self.logfile and not self.closed:
            self.logfile.close()
            self.closed = True

    def info(self, *args, **kwargs):
        self.log(*args, level=0, **kwargs)

    def warn(self, *args, **kwargs):
        self.log(*args, level=1, **kwargs)

    def error(self, *args, **kwargs):
        self.log(*args, level=2, **kwargs)

    def log(self, *args, level=0, **kwargs):
        """
        Log output to STDERR and/or the specified log file
        :param args: arguments to log
        :param level: 0 = INFO, 1 = WARNING, 2 = ERROR
        :param kwargs: kwargs to send to the print() statement
        :return:
        """

        if level == 0:
            self.infos += 1
        elif level == 1:
            self.warnings += 1
        elif level == 2:
            self.errors += 1

        if level >= self.level:
            if level == 0:
                key = 'INFO:'
            elif level == 1:
                key = 'WARN:'
            elif level == 2:
                key = 'ERROR:'
            else:
                key = None
            if self.console:
                print(key, *args, file=self.console, **kwargs)
            if self.logfile:
                print(key, *args, file=self.logfile, **kwargs)

    def get_num_messages(self):
        return self.infos, self.warnings, self.errors


def gen_sha256(file):
    sha256_hash = hashlib.sha256()
    with open(file, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def to_datetime(datetime_str):
    datetime_object = datetime.datetime.strptime(datetime_str, _DATETIME_FORMAT)
    return datetime_object


def read_config_file(filename):

    if os.path.isfile(filename):
        if filename.endswith('.yaml'):
            with open(filename, 'r') as stream:
                config = yaml.safe_load(stream)
                return config
        elif filename.endswith('.json'):
            with open(filename, 'r') as stream:
                config = json.load(stream)
                return config
        else:
            raise ValueError('Only .json or .yaml files are supported. {} was specified'.format(filename))
    else:
        msg = 'Config file {} not found'.format(filename)
        raise ValueError(msg)

def find_property(my_dict, key, default=None):
    if key in my_dict:
        return my_dict[key]
    else:
        return default


def find_path(my_dict, key, default=None):
    value = find_property(my_dict, key, default=default)
    if value:
        return Path(value)
    else:
        return default


def make_path_relative(p):
    if p.is_absolute():
        return p.relative_to('/')
    else:
        return p


def expand_path(p1, p2, expand=True):
    if expand and p1:
        return p1 / make_path_relative(p2)
    else:
        return p2


def main():

    log = Logger(logfile='logfile.log', level=1)

    log.log('a', 'b', 'c', level=0)
    log.log('foo', 'bar', 'baz')
    log.log('foo', 99, 'apples', level=2)


if __name__ == "__main__":
    main()
