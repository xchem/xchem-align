import os
from pathlib import Path
import shutil

import pytest

@pytest.fixture(scope="session")
def constants():
    class Constants:
        TEST_DATA_DIR = "test-data"
        INPUT_1_DIR = "test-data/input_1"
        INPUT_2_DIR = "test-data/input_2"
        CONFIG_1_FILE = "test-data/config_1.yaml"
        CONFIG_2_FILE = "test-data/config_2.yaml"
        ASSEMBLIES_FILE = "test-data/outputs/assemblies.yaml"
        XTAFLORMS_FILE = "test-data/outputs/crystalforms.yaml"
        TEST_DIR = "test-data/outputs"
        UPLOAD_1_DIR = "test-data/outputs/upload_1"
        UPLOAD_2_DIR = "test-data/outputs/upload_2"
        METADATA_FILE = "meta_collator.yaml"

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

    path = Path(constants.TEST_DIR)
    return path


@pytest.fixture(scope="session")
def upload_1_dir(constants, test_dir):
    path = Path(constants.UPLOAD_1_DIR)
    if path.exists():
        shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def upload_2_dir(constants, test_dir):
    path = Path(constants.UPLOAD_2_DIR)
    if path.exists():
        shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def upload_1_data_dir(constants, test_dir):
    path = Path(constants.INPUT_1_DIR)
    return path


@pytest.fixture(scope="session")
def upload_2_data_dir(constants, test_dir):
    path = Path(constants.INPUT_2_DIR)
    return path

@pytest.fixture(scope="session")
def config_1_file(constants, ):
    path = Path(constants.CONFIG_1_FILE)
    return path

@pytest.fixture(scope="session")
def config_2_file(constants, ):
    path = Path(constants.CONFIG_2_FILE)
    return path

@pytest.fixture(scope="session")
def assemblies_file(constants, ):
    path = Path(constants.ASSEMBLIES_FILE)
    return path

@pytest.fixture(scope="session")
def xtalforms_file(constants, ):
    path = Path(constants.XTAFLORMS_FILE)
    return path
