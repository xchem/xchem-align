import os
from pathlib import Path
import shutil

import yaml
import pytest

from ligand_neighbourhood_alignment import dt


@pytest.fixture(scope="session")
def constants():
    class Constants:
        TEST_DATA_DIR = "test-data"
        INPUT_1_DIR = "test-data/input_1"
        INPUT_2_DIR = "test-data/input_2"
        INPUT_3_DIR = "test-data/input_3"
        CONFIG_1_FILE = "test-data/config_1.yaml"
        CONFIG_2_FILE = "test-data/config_2.yaml"
        CONFIG_3_FILE = "test-data/config_3.yaml"
        TEST_DIR = "test-data/outputs"
        CONFIG_FILE = TEST_DIR + '/upload-current/config.yaml'
        UPLOAD_1_DIR = TEST_DIR + '/upload-current/upload_1'
        UPLOAD_2_DIR = TEST_DIR + '/upload-current/upload_2'
        UPLOAD_3_DIR = TEST_DIR + '/upload-current/upload_3'
        METADATA_FILE = "meta_collator.yaml"

        # TEST_DATA_DIR = "data"
        ASSEMBLIES_FILE = "test-data/assemblies.yaml"
        TEST_OUTPUT_DIR = "tests/output"
        DATA_PATHS = {
            "5rgs": "test-data/5rgs.pdb",
            "8e1y": "test-data/8e1y.pdb",
            "8dz9": "test-data/8dz9.pdb",
            "7ql8": "test-data/7ql8-pdb-bundle1.pdb",
            "Mpro-i0130": "test-data/Mpro-i0130.pdb",
            "Mpro-IBM0078": "test-data/refine_6.split.bound-state.pdb",
            "Mpro-IBM0058": "test-data/refine_7.split.bound-state.pdb",
            "Mpro-x0107": "test-data/refine_8.split.bound-state.pdb",
            "Mpro-IBM0045": "test-data/refine_16.split.bound-state.pdb",
        }

    return Constants()


@pytest.fixture(scope="session")
def test_data_dir(constants):
    return Path(constants.TEST_DATA_DIR)


@pytest.fixture(scope="session")
def test_dir(constants):
    path = Path(constants.UPLOAD_1_DIR)
    if path.exists():
        shutil.rmtree(path)

    path = Path(constants.UPLOAD_2_DIR)
    if path.exists():
        shutil.rmtree(path)

    path = Path(constants.UPLOAD_3_DIR)
    if path.exists():
        shutil.rmtree(path)

    path = Path(constants.TEST_DIR)
    return path


@pytest.fixture(scope="session")
def upload_1_dir(constants, test_dir):
    for path in [Path(constants.UPLOAD_3_DIR), Path(constants.UPLOAD_2_DIR), Path(constants.UPLOAD_1_DIR)]:
        if path.exists():
            shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def upload_2_dir(constants, test_dir):
    for path in [Path(constants.UPLOAD_3_DIR), Path(constants.UPLOAD_2_DIR)]:
        if path.exists():
            shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def upload_3_dir(constants, test_dir):
    path = Path(constants.UPLOAD_3_DIR)
    if path.exists():
        shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def config_1_file(
    constants,
):
    path = Path(constants.CONFIG_FILE)
    shutil.copy(constants.CONFIG_1_FILE, constants.CONFIG_FILE)
    return path


@pytest.fixture(scope="session")
def config_2_file(
    constants,
):
    path = Path(constants.CONFIG_FILE)
    shutil.copy(constants.CONFIG_2_FILE, constants.CONFIG_FILE)
    return path


@pytest.fixture(scope="session")
def config_3_file(
    constants,
):
    path = Path(constants.CONFIG_FILE)
    shutil.copy(constants.CONFIG_3_FILE, constants.CONFIG_FILE)
    return path


########################
### from lna package ###
########################


# @pytest.fixture(scope="session")
# def test_data_dir(constants):
#     return Path(constants.TEST_DATA_DIR)


@pytest.fixture(scope="session")
def test_output_dir(constants):
    path = Path(constants.TEST_OUTPUT_DIR)
    if path.exists():
        shutil.rmtree(path)

    return path


@pytest.fixture(scope="session")
def assemblies_file(
    constants,
):
    path = Path(constants.ASSEMBLIES_FILE)
    return path


@pytest.fixture(scope="session")
def pdb_paths(constants):
    pdb_paths = {key: Path(path) for key, path in constants.DATA_PATHS.items()}
    return pdb_paths


@pytest.fixture(scope="session")
def assemblies(constants, assemblies_file):
    _assemblies = {}
    with open(assemblies_file, "r") as f:
        dic = yaml.safe_load(f)

    for assembly_id, assembly_info in dic["assemblies"].items():
        _assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return _assemblies


@pytest.fixture(scope="session")
def xtalforms(constants, assemblies_file):
    _xtalforms = {}
    with open(assemblies_file, "r") as f:
        dic = yaml.safe_load(f)["crystalforms"]

    for xtalform_id, xtalform_info in dic.items():
        _xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return _xtalforms
