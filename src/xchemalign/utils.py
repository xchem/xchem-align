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

import atexit
import datetime
import hashlib
import os
from pathlib import Path
import sys
import json

import yaml

from rdkit import Chem, Geometry
from gemmi import cif

_DATETIME_FORMAT = "%Y-%m-%d %H:%M:%S"


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
        "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.native.ccp4",
    ]
    ASSEMBLIES_FILENAME = "assemblies.yaml"
    METADATA_XTAL_FILENAME = "meta_collator{}.yaml"
    METADATA_ALIGN_FILENAME = "meta_aligner.yaml"
    VERSION_DIR_PREFIX = "upload_"
    DEFAULT_SOAKDB_PATH = "processing/database/soakDBDataFile.sqlite"
    DEFAULT_MODEL_BUILDING_DIR = "processing/analysis/model_building"
    CONFIG_INPUTS = "inputs"
    CONFIG_TYPE = "type"
    CONFIG_DIR = "dir"
    CONFIG_CWD = "collator_cwd"
    CONFIG_CONFIG_FILE = "collator_config_file"
    CONFIG_SOAKDB = "soakdb"
    CONFIG_TYPE_MODEL_BUILDING = "model_building"
    CONFIG_TYPE_MANUAL = "manual"
    CONFIG_BASE_DIR = "base_dir"
    CONFIG_OUTPUT_DIR = "output_dir"
    CONFIG_TARGET_NAME = "target_name"
    CONFIG_REF_DATASETS = "ref_datasets"
    CONFIG_EXCLUDE = 'exclude'
    CONFIG_CODE_PREFIX = 'code_prefix'
    CONFIG_CODE_PREFIX_TOOLTIP = 'code_prefix_tooltip'
    META_RUN_ON = "run_on"
    META_INPUT_DIRS = "input_dirs"
    META_VERSION_NUM = "version_number"
    META_VERSION_DIR = "version_dir"
    META_PREV_VERSION_DIRS = "previous_version_dirs"
    META_LAST_UPDATED = "last_updated"
    META_REFINEMENT_OUTCOME = "refinement_outcome"
    META_STATUS = "status"
    META_STATUS_SUPERSEDED = "superseded"
    META_STATUS_UNCHANGED = "unchanged"
    META_STATUS_NEW = "new"
    META_STATUS_DEPRECATED = "deprecated"
    META_STATUS_REASON = "status_reason"
    META_XTALS = "crystals"
    META_REFERENCE = "reference"
    META_FILE = "file"
    META_SHA256 = "sha256"
    META_XTAL_FILES = "crystallographic_files"
    META_ALIGNED_FILES = "aligned_files"
    META_REFERENCE_ALIGNMENTS = "reference_aligned_files"
    META_XTAL_PDB = "xtal_pdb"
    META_XTAL_MTZ = "xtal_mtz"
    META_XTAL_CIF = "ligand_cif"
    META_SMILES = "smiles"
    META_BINDING_EVENT = "panddas_event_files"
    META_PANDDAS_MISSING_OK = "panddas_missing_ok"
    META_PROT_MODEL = "model"
    META_PROT_CHAIN = "chain"
    META_PROT_RES = "res"
    META_PROT_INDEX = "index"
    META_PROT_BDC = "bdc"
    META_AIGNED_STRUCTURE = "structure"
    META_AIGNED_ARTEFACTS = "artefacts"
    META_AIGNED_EVENT_MAP = "event_map"
    META_AIGNED_X_MAP = "sigmaa_map"
    META_AIGNED_DIFF_MAP = "diff_map"
    META_CONFORMER_SITES = "conformer_sites"
    META_CONFORMER_SITE_NAME = "name"
    META_CONFORMER_SITE_REFERENCE_LIG = "lig_ref"
    META_CONFORMER_SITE_RESIDUES = "residues"
    META_CONFORMER_SITE_MEMBERS = "members"
    META_CANONICAL_SITES = "canon_sites"
    META_CANONICAL_SITE_REF_SUBSITE = "canon_site_ref_site"
    META_CANONICAL_SITE_CONFORMER_SITES = "canon_site_conf_sites"
    META_CANONICAL_SITE_RESIDUES = "site_residues"
    META_CANONICAL_SITE_MEMBERS = "site_members"
    META_XTALFORM_SITES = "xtalform_sites"
    META_XTALFORM_SITE_XTALFORM_ID = "xtalform_id"
    META_XTALFORM_SITE_CANONICAL_SITE_ID = "canon_site_id"
    META_XTALFORM_SITE_LIGAND_CHAIN = "lig_chain"
    META_XTALFORM_SITE_MEMBERS = "xtalform_members"
    META_ASSIGNED_XTALFORM = "assigned_xtalform"
    META_CHAIN = "chain"
    META_RESIDUE = "residue"
    META_DTAG = "dtag"
    META_XTALFORM_REFERENCE = "xtalform_ref"
    META_XTALFORM_SPACEGROUP = "xtalform_space_group"
    META_XTALFORM_CELL = "xtalform_cell"
    META_ASSEMBLIES_XTALFORMS = "assemblies_xtalforms"
    META_XTALFORMS = "crystalforms"
    META_ASSEMBLIES = "assemblies"
    META_ASSEMBLIES_CHAINS = "chains"
    META_PDB_APO = "pdb_apo"
    META_PDB_APO_SOLV = "pdb_apo_solv"
    META_PDB_APO_DESOLV = "pdb_apo_desolv"
    META_LIGAND_MOL = "ligand_mol"
    META_LIGAND_PDB = "ligand_pdb"
    META_LIGAND_SMILES_STRING = "ligand_smiles_string"
    META_LIGAND_SMILES = "ligand_smiles"
    META_TRANSFORMS = "transforms"
    META_TRANSFORMS_OBSERVATION_TO_CONFORMER_SITES = "observation_to_conformer"
    META_TRANSFORMS_CONFORMER_SITES_TO_CANON = "conformer_to_canon"
    META_TRANSFORMS_CANON_SITES_TO_GLOBAL = "canon_to_global"
    META_TRANSFORMS_GLOBAL_REFERENCE_CANON_SITE_ID = "global_reference_canon_site_id"
    META_CMPD_CODE = "compound_code"
    META_GIT_INFO = "xca_git_info"
    META_GIT_INFO_URL = "origin_url"
    META_GIT_INFO_BRANCH = "branch"
    META_GIT_INFO_SHA = "sha"
    META_GIT_INFO_TAG = "tag"
    META_GIT_INFO_DIRTY = "dirty"
    META_CODE_PREFIX = "code_prefix"
    META_CODE_PREFIX_TOOLTIPS = "code_prefix_tooltips"
    SOAKDB_XTAL_NAME = "CrystalName"
    SOAKDB_COL_PDB = "RefinementBoundConformation"
    SOAKDB_COL_MTZ = "RefinementMTZ_latest"
    SOAKDB_COL_CIF = "RefinementCIF"
    SOAKDB_COL_LAST_UPDATED = "LastUpdatedDate"
    SOAKDB_COL_REFINEMENT_OUTCOME = "RefinementOutcome"
    SOAKDB_COL_COMPOUND_CODE = "CompoundCode"
    CRYSTAL_NEW = "crystal_new"
    ASSEMBLIES_FILENAME = "assemblies.yaml"
    PREVIOUS_OUTPUT_DIR = ""
    ENV_XCA_GIT_REPO = "XCA_GIT_REPO"


BOND_TYPES = {
    'single': Chem.rdchem.BondType.SINGLE,
    'double': Chem.rdchem.BondType.DOUBLE,
    'triple': Chem.rdchem.BondType.TRIPLE,
    'SINGLE': Chem.rdchem.BondType.SINGLE,
    'DOUBLE': Chem.rdchem.BondType.DOUBLE,
    'TRIPLE': Chem.rdchem.BondType.TRIPLE,
    'aromatic': Chem.rdchem.BondType.AROMATIC,
    'deloc': Chem.rdchem.BondType.SINGLE,
}


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
        self.infos = []
        self.warnings = []
        self.errors = []
        if logfile:
            self.logfile = open(logfile, "w")
            self.logfilename = logfile
            self.closed = False
        else:
            self.logfile = None
            self.logfilename = None
            self.closed = True
        atexit.register(self.close)
        x = datetime.datetime.now()
        self.log("initialising logging at level {} at {}".format(level, x), level=0)
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

        msg = " ".join([str(s) for s in args])

        if level == 0:
            self.infos.append(msg)
        elif level == 1:
            self.warnings.append(msg)
        elif level == 2:
            self.errors.append(msg)

        if level >= self.level:
            if level == 0:
                key = "INFO:"
            elif level == 1:
                key = "WARN:"
            elif level == 2:
                key = "ERROR:"
            else:
                key = None
            if self.console:
                print(key, *args, file=self.console, **kwargs)
            if self.logfile:
                print(key, *args, file=self.logfile, **kwargs)

    def get_num_messages(self):
        return len(self.infos), len(self.warnings), len(self.errors)

    def report(self):
        """
        Write out a summary of the warnings and errors
        :return:
        """
        if self.warnings:
            print(
                "*********** CAUTION:",
                len(self.warnings),
                "warnings were generated. See above for context ***********",
            )
            for msg in self.warnings:
                print("WARN:", msg)

        if self.errors:
            print("*********** CAUTION:", len(self.errors), "errors were generated. See above for context ***********")
            for msg in self.errors:
                print("ERROR:", msg)


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
        if filename.endswith(".yaml"):
            with open(filename, "r") as stream:
                config = yaml.safe_load(stream)
                return config
        elif filename.endswith(".json"):
            with open(filename, "r") as stream:
                config = json.load(stream)
                return config
        else:
            raise ValueError("Only .json or .yaml files are supported. {} was specified".format(filename))
    else:
        msg = "Config file {} not found".format(filename)
        raise ValueError(msg)


def find_property(my_dict, key, default=None):
    if key in my_dict:
        v = my_dict[key]
        # value can be None if the YAML tag is defined but has no values
        if v is None:
            return default
        else:
            return v
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
        return p.relative_to("/")
    else:
        return p


def expand_path(p1, p2, expand=True):
    if expand and p1:
        return p1 / make_path_relative(p2)
    else:
        return p2


def gen_mol_from_cif(cif_file):
    mol = Chem.RWMol()
    conf = Chem.Conformer()

    doc = cif.read(cif_file)
    # Diamond CIFs have two blocks, but the one we want will be named data_comp_LIG
    block = doc.find_block("comp_LIG")
    # Other CIFs have unpredictable block names, so let's hope there is only one
    if not block:
        block = doc.sole_block()
    if not block:
        print("sole block not found")
        return None
    atom_ids = block.find_loop('_chem_comp_atom.atom_id')
    atom_symbols = block.find_loop('_chem_comp_atom.type_symbol')
    # coordinates are sometimes called "x" and sometimes "model_Cartn_x" etc.
    x = block.find_loop('_chem_comp_atom.x')
    if not x:
        x = block.find_loop('_chem_comp_atom.model_Cartn_x')
    y = block.find_loop('_chem_comp_atom.y')
    if not y:
        y = block.find_loop('_chem_comp_atom.model_Cartn_y')
    z = block.find_loop('_chem_comp_atom.z')
    if not z:
        z = block.find_loop('_chem_comp_atom.model_Cartn_z')
    charges = [0] * len(atom_ids)
    if block.find_loop('_chem_comp_atom.charge'):
        charges = list(block.find_loop('_chem_comp_atom.charge'))
    elif block.find_loop('_chem_comp_atom.partial_charge'):
        charges = list(block.find_loop('_chem_comp_atom.partial_charge'))

    atoms = {}
    for s, id, px, py, pz, charge in zip(atom_symbols, atom_ids, x, y, z, charges):
        # sometimes that atom ids are wrapped in double quotes
        if id[0] == '"':
            id = id[1:]
        if id[-1] == '"':
            id = id[:-1]

        if len(s) == 2:
            s = s[0] + s[1].lower()

        atom = Chem.Atom(s)
        atom.SetFormalCharge(round(float(charge)))
        atom.SetProp('atom_id', id)
        idx = mol.AddAtom(atom)
        atom.SetIntProp('idx', idx)
        atoms[id] = atom

        point = Geometry.Point3D(float(px), float(py), float(pz))
        conf.SetAtomPosition(idx, point)

    atom1 = block.find_loop('_chem_comp_bond.atom_id_1')
    atom2 = block.find_loop('_chem_comp_bond.atom_id_2')
    bond_type = block.find_loop('_chem_comp_bond.type')

    for a1, a2, bt in zip(atom1, atom2, bond_type):
        mol.AddBond(atoms[a1].GetIntProp('idx'), atoms[a2].GetIntProp('idx'), BOND_TYPES[bt])

    Chem.SanitizeMol(mol)
    mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol)
    mol = Chem.RemoveAllHs(mol)

    return mol


def main():
    # log = Logger(logfile="logfile.log", level=1)
    #
    # log.log("a", "b", "c", level=0)
    # log.log("foo", "bar", "baz")
    # log.log("foo", 99, "apples", level=2)

    mol = gen_mol_from_cif('data/Zx1674a.cif')
    molfile = Chem.MolToMolBlock(mol)

    print(molfile)


if __name__ == "__main__":
    main()
