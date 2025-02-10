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
