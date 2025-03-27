from pathlib import Path

import pytest
import yaml

from xchemalign.collator import Collator
from xchemalign.aligner import Aligner
from xchemalign import utils
from xchemalign.utils import Constants


def test_collator_upload_1(constants, test_dir, uploads_dir, config_1_file, upload_1_dir, uploads_dir):
    c = Collator(uploads_dir)
    logger = c.logger

    meta = c.validate()

    if meta is None or len(logger.errors) > 0:
        print("There are errors, cannot continue")
        print(logger.errors)
        exit(1)
    else:
        c.run(meta)


@pytest.mark.order(after="test_collator_upload_1")
def test_aligner_upload_1(constants):
    log = str(Path(constants.TEST_DIR) / 'aligner.log')

    a = Aligner(constants.TEST_DIR, log_file=log, log_level=0)
    logger = a.logger
    utils.LOG = logger

    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        print(a.logger.errors)
        exit(1)
    else:
        a.run()


@pytest.mark.order(after="test_aligner_upload_1")
def test_collator_upload_2(constants, config_2_file, upload_2_dir, uploads_dir):
    c = Collator(uploads_dir)
    logger = c.logger

    meta = c.validate()

    print(meta)

    if meta is None or len(logger.errors) > 0:
        print("There are errors, cannot continue")
        print(logger.errors)
        exit(1)
    else:
        c.run(meta)
    assert (
        # len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {}))
        len(
            meta[Constants.META_XTALS]["Mpro-x0107_fake_P1"][Constants.META_XTAL_FILES].get(
                Constants.META_BINDING_EVENT, {}
            )
        )
        != 0
    )


@pytest.mark.order(after="test_collator_upload_2")
def test_aligner_upload_2(constants):
    log = str(Path(constants.TEST_DIR) / 'aligner.log')

    a = Aligner(constants.TEST_DIR, log_file=log, log_level=0)
    logger = a.logger
    utils.LOG = logger

    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        print(a.logger.errors)
        exit(1)
    else:
        a.run()
    assert "Mpro-i0130" in [x.name for x in (Path(constants.UPLOAD_2_DIR) / "aligned_files").glob("*")]


# @pytest.mark.order(after="test_aligner_upload_2")
# def test_collator_upload_3(
#     constants,
#     test_dir,
#     upload_3_dir,
#     upload_3_data_dir,
#     config_3_file,
# ):
#     c = Collator(str(config_3_file), )
#
#     meta, num_errors, num_warnings = c.validate()
#
#     if meta is None or num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         c.run(meta)
#
#     with open(Path(upload_3_dir) / "meta_collator.yaml", 'r') as f:
#         new_meta = yaml.safe_load(f)
#     assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "crystallographic_files").glob("*")]
#
#
# @pytest.mark.order(after="test_collator_upload_3")
# def test_aligner_upload_3(constants, xtalforms_file, upload_3_dir):
#     a = Aligner(upload_3_dir, constants.METADATA_FILE, xtalforms_file, )
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         a.run()
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "aligned_files").glob("*")]
