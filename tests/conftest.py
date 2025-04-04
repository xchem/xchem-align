import os
import subprocess
from pathlib import Path
import shutil

import pytest


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
        (Path(constants.TEST_DIR) / constants.VERSION_DIR ).resolve(),
        (Path(constants.TEST_DIR) / 'upload-current').resolve(),
        target_is_directory=True)

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
    path = Path(Path(constants.TEST_DIR) / constants.VERSION_DIR  / constants.CONFIG_YAML)
    shutil.copy(Path(constants.TEST_DATA_DIR) / constants.CONFIG_1_FILE, path)
    return path


@pytest.fixture(scope="session")
def config_2_file(
    constants,
    uploads_dir
):
    # os.remove(Path(constants.TEST_DIR) / 'upload-current')
    # os.symlink(
    #     (Path(constants.TEST_DIR) / constants.VERSION_DIR / 'upload_2').resolve(),
    #     (Path(constants.TEST_DIR) / 'upload-current').resolve(),
    #
    #     target_is_directory=True)
    path = Path(Path(constants.TEST_DIR) / constants.VERSION_DIR  / constants.CONFIG_YAML)
    os.remove(path)
    shutil.copy(Path(constants.TEST_DATA_DIR) / constants.CONFIG_2_FILE, path)
    return path


# @pytest.fixture(scope="session")
# def config_3_file(
#     constants,
# ):
#     path = Path(constants.CONFIG_FILE)
#     shutil.copy(constants.CONFIG_3_FILE, constants.CONFIG_FILE)
#     return path

