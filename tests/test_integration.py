from pathlib import Path

import pytest
import yaml

from xchemalign.collator import Collator
from xchemalign.aligner import Aligner
from xchemalign.utils import Constants

def test_collator_upload_1(
        constants,
        test_dir,
        upload_1_dir,
        upload_1_data_dir,
        config_1_file,
):
    c = Collator(str(config_1_file), )
    logger = c.logger

    meta = c.validate()

    if meta is None or len(logger.errors)>0:
        print("There are errors, cannot continue")
        exit(1)
    else:
        c.run(meta)



@pytest.mark.order(after="test_collator_upload_1")
def test_aligner_upload_1(constants, xtalforms_file, upload_1_dir):
    a = Aligner(upload_1_dir, constants.METADATA_FILE, xtalforms_file, )
    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        a.run()



@pytest.mark.order(after="test_aligner_upload_1")
def test_collator_upload_2(
    constants,
    test_dir,
    upload_2_dir,
    upload_2_data_dir,
    config_2_file,
):
    c = Collator(str(config_2_file), )
    logger = c.logger

    meta = c.validate()

    if meta is None or len(logger.errors)>0:
        print("There are errors, cannot continue")
        exit(1)
    else:
        c.run(meta)

    with open(Path(upload_2_dir) / "meta_collator.yaml", 'r') as f:
        new_meta = yaml.safe_load(f)
    assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0


@pytest.mark.order(after="test_collator_upload_2")
def test_aligner_upload_2(constants, xtalforms_file, upload_2_dir):
    a = Aligner(upload_2_dir, constants.METADATA_FILE, xtalforms_file, )
    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        a.run()

    assert "Mpro-i0130" in [x.name for x in (Path(upload_2_dir) / "aligned_files").glob("*")]

@pytest.mark.order(after="test_aligner_upload_2")
def test_collator_upload_2(
    constants,
    test_dir,
    upload_3_dir,
    upload_3_data_dir,
    config_3_file,
):
    c = Collator(str(config_3_file), )

    meta, num_errors, num_warnings = c.validate()

    if meta is None or num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        c.run(meta)

    with open(Path(upload_3_dir) / "meta_collator.yaml", 'r') as f:
        new_meta = yaml.safe_load(f)
    assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0

    assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "crystallographic_files").glob("*")]


@pytest.mark.order(after="test_collator_upload_3")
def test_aligner_upload_3(constants, xtalforms_file, upload_3_dir):
    a = Aligner(upload_3_dir, constants.METADATA_FILE, xtalforms_file, )
    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        a.run()

    assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "aligned_files").glob("*")]

