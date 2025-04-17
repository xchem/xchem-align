import os
import subprocess
from pathlib import Path
import shutil

import yaml
import pytest

from ligand_neighbourhood_alignment import dt


@pytest.fixture(scope="session")
def constants():
    class Constants:
        TEST_DATA_DIR = "test-data"
        VERSION_DIR = "upload-v2"
        INPUT_1_DIR = "test-data/input_1"
        INPUT_2_DIR = "test-data/input_2"
        INPUT_3_DIR = "test-data/input_3"
        CONFIG_YAML = "config.yaml"
        CONFIG_1_FILE = "config_1.yaml"
        CONFIG_2_FILE = "config_2.yaml"
        CONFIG_3_FILE = "test-data/config_3.yaml"
        TEST_DIR = "test-data/outputs"
        CONFIG_FILE = TEST_DIR + '/upload-v2/config.yaml'
        UPLOAD_1_DIR = TEST_DIR + '/upload-v2/upload_1'
        UPLOAD_2_DIR = TEST_DIR + '/upload-v2/upload_2'
        UPLOAD_3_DIR = TEST_DIR + '/upload-v2/upload_3'
        METADATA_FILE = "meta_collator.yaml"
        ASSEMBLIES_FILENAME = "assemblies.yaml"

        # TEST_DATA_DIR = "data"
        ASSEMBLIES_FILE = "test-data/assemblies.yaml"

        LNA_ASSEMBLIES_FILENAME = "lna_assemblies.yaml"
        LNA_ASSEMBLIES_FILE = "test-data/lna_assemblies.yaml"

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
    path = Path(constants.TEST_DIR) / constants.VERSION_DIR
    if path.exists():
        shutil.rmtree(path)

    path = (Path(constants.TEST_DIR) / 'upload-current').resolve()
    if path.exists():
        shutil.rmtree(path)
        os.unlink(path)

    # path = Path(constants.UPLOAD_2_DIR)
    # if path.exists():
    #     shutil.rmtree(path)
    #
    # path = Path(constants.UPLOAD_3_DIR)
    # if path.exists():
    #     shutil.rmtree(path)

    path = Path(constants.TEST_DIR)
    return path


@pytest.fixture(scope="session")
def uploads_dir(constants, test_dir):
    os.mkdir(Path(constants.TEST_DIR) / constants.VERSION_DIR)
    os.mkdir(Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_1')
    # os.mkdir(Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload-current')

    if (Path(constants.TEST_DIR) / 'upload-current').resolve().exists():
        # p = subprocess.Popen(f"unlink {(Path(constants.TEST_DIR) / 'upload-current')}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # stdout, stderr = p.communicate()
        os.remove((Path(constants.TEST_DIR) / 'upload-current'))
        # print(stdout)
        # print(stderr)

    os.symlink(
        (Path(constants.TEST_DIR) / constants.VERSION_DIR).resolve(),
        (Path(constants.TEST_DIR) / 'upload-current').resolve(),
        target_is_directory=True,
    )

    return Path(constants.TEST_DIR) / constants.VERSION_DIR


@pytest.fixture(scope="session")
def upload_1_dir(constants, test_dir, uploads_dir):
    # for path in [Path(constants.UPLOAD_3_DIR), Path(constants.UPLOAD_2_DIR), Path(constants.UPLOAD_1_DIR)]:
    #     if path.exists():
    #         shutil.rmtree(path)
    # os.mkdir(path)
    path = Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_1'
    return path


@pytest.fixture(scope="session")
def upload_2_dir(constants, test_dir, uploads_dir):
    # for path in [Path(constants.UPLOAD_3_DIR), Path(constants.UPLOAD_2_DIR)]:
    #     if path.exists():
    #         shutil.rmtree(path)
    # os.mkdir(path)
    os.mkdir(Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_2')
    path = Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_2'

    return path


# @pytest.fixture(scope="session")
# def upload_3_dir(constants, test_dir):
#     path = Path(constants.UPLOAD_3_DIR)
#     if path.exists():
#         shutil.rmtree(path)
#     os.mkdir(path)
#
#     return path


@pytest.fixture(scope="session")
def assemblies_file(
    constants,
    uploads_dir,
):
    path = Path(Path(constants.TEST_DIR) / constants.VERSION_DIR / constants.ASSEMBLIES_FILENAME)
    shutil.copy(Path(constants.TEST_DATA_DIR) / constants.ASSEMBLIES_FILENAME, path)
    return path


@pytest.fixture(scope="session")
def config_1_file(
    constants,
    uploads_dir,
):
    path = Path(Path(constants.TEST_DIR) / constants.VERSION_DIR / constants.CONFIG_YAML)
    shutil.copy(Path(constants.TEST_DATA_DIR) / constants.CONFIG_1_FILE, path)
    return path


@pytest.fixture(scope="session")
def config_2_file(constants, uploads_dir):
    # os.remove(Path(constants.TEST_DIR) / 'upload-current')
    # os.symlink(
    #     (Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_2').resolve(),
    #     (Path(constants.TEST_DIR) / 'upload-current').resolve(),
    #
    #     target_is_directory=True)
    path = Path(Path(constants.TEST_DIR) / constants.VERSION_DIR / constants.CONFIG_YAML)
    os.remove(path)
    shutil.copy(Path(constants.TEST_DATA_DIR) / constants.CONFIG_2_FILE, path)
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


# TODO: probably just needs a rename
# @pytest.fixture(scope="session")
# def assemblies_file(
#     constants,
# ):
#     path = Path(constants.ASSEMBLIES_FILE)
#     return path


# rename attempt
@pytest.fixture(scope="session")
def lna_assemblies_file(
    constants,
):
    path = Path(constants.LNA_ASSEMBLIES_FILE)
    return path


@pytest.fixture(scope="session")
def pdb_paths(constants):
    pdb_paths = {key: Path(path) for key, path in constants.DATA_PATHS.items()}
    return pdb_paths


@pytest.fixture(scope="session")
def lna_assemblies(constants, lna_assemblies_file):
    _assemblies = {}
    with open(lna_assemblies_file, "r") as f:
        dic = yaml.safe_load(f)

    for assembly_id, assembly_info in dic["assemblies"].items():
        _assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return _assemblies


@pytest.fixture(scope="session")
def xtalforms(constants, lna_assemblies_file):
    _xtalforms = {}
    with open(lna_assemblies_file, "r") as f:
        dic = yaml.safe_load(f)["crystalforms"]

    for xtalform_id, xtalform_info in dic.items():
        _xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return _xtalforms
