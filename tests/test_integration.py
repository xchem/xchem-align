from pathlib import Path

import pytest

from xchemalign.collator import Collator
from xchemalign.aligner import Aligner

def test_collator_upload_1(
        constants,
        test_dir,
        upload_1_dir,
        upload_1_data_dir,
        config_1_file,
):
    c = Collator(str(config_1_file), )

    meta, num_errors, num_warnings = c.validate()

    if meta is None or num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        c.run(meta)



@pytest.mark.order(after="test_collator_upload_1")
def test_aligner_upload_1(constants, assemblies_file, xtalforms_file, upload_1_dir):
    a = Aligner(upload_1_dir, constants.METADATA_FILE, xtalforms_file, assemblies_file, )
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

    meta, num_errors, num_warnings = c.validate()

    if meta is None or num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        c.run(meta)

@pytest.mark.order(after="test_collator_upload_2")
def test_aligner_upload_2(constants, assemblies_file, xtalforms_file, upload_2_dir):
    a = Aligner(upload_2_dir, constants.METADATA_FILE, xtalforms_file, assemblies_file, )
    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        exit(1)
    else:
        a.run()

    assert "Mpro-i0130" in [x.name for x in (Path(upload_2_dir) / "aligned_files").glob("*")]
