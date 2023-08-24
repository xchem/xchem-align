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
        CONFIG_1_FILE = "test-data/config.yaml"
        CONFIG_2_FILE = "test-data/config_2.yaml"
        ASSEMBLIES_FILE = "test-data/outputs/assemblies.yaml"
        XTAFLORMS_FILE = "test-data/outputs/xtalforms.yaml"
        TEST_DIR = "test-data/outputs"
        UPLOAD_1_DIR = "test-data/outputs/upload_1"
        UPLOAD_2_DIR = "test-data/outputs/upload_2"

    return Constants()

@pytest.fixture(scope="session")
def test_data_dir(constants):
@pytest.fixture(scope="session")
def test_dir(constants):
    path = Path(constants.TEST_DIR)
    if path.exists():
        shutil.rmtree(path)
    os.mkdir(path)

    return path


@pytest.fixture(scope="session")
def upload_1_dir(constants, test_dir):
    ...


@pytest.fixture(scope="session")
def upload_2_dir(constants, test_dir):
    ...


@pytest.fixture(scope="session")
def upload_1_data_dir(constants, test_dir):
    ...


@pytest.fixture(scope="session")
def upload_2_data_dir(constants, test_dir):
    ...
