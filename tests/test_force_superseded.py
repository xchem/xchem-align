"""
Integration test for the collator's ``--force-superseded`` option (issue #87).

Scenario:
  * upload_1 is run with one dataset (Mpro-IBM0078) excluded.
  * upload_2 is run with that dataset included and with ``force_superseded=True``.

Expected outcome:
  * the previously excluded dataset appears in upload_2 with status ``new``;
  * datasets present in both uploads (unchanged PDB) are forced to status ``superseded``;
  * the forced datasets have their crystallographic files re-copied into the upload_2
    directory and their metadata file paths point at upload_2 (not the previous upload).
"""

import os
import shutil
from pathlib import Path

import yaml

from xchemalign.collator import Collator
from xchemalign.utils import Constants


# dataset excluded from upload_1 then re-introduced in upload_2 -> should become "new"
EXCLUDED_DATASET = "Mpro-IBM0078"
# datasets present (unchanged) in both uploads -> forced to "superseded" in upload_2
UNCHANGED_DATASETS = ["Mpro-IBM0045", "Mpro-IBM0058"]


def _read_meta(upload_dir: Path) -> dict:
    with open(upload_dir / Constants.METADATA_XTAL_FILENAME.format("")) as f:
        return yaml.safe_load(f)


def _setup_working_dir(tmp_path: Path, exclude: list[str]) -> tuple[Path, Path]:
    """Create a working dir with an upload-current symlink, assemblies.yaml and a config.yaml.

    The config is based on test-data/config_1.yaml; ``exclude`` is applied to the
    model_building input. Returns (working_dir, version_dir).
    """
    version_dir = tmp_path / "upload-v3"
    (version_dir / "upload_1").mkdir(parents=True)
    os.symlink(version_dir, tmp_path / "upload-current", target_is_directory=True)

    shutil.copy("test-data/assemblies.yaml", version_dir / "assemblies.yaml")

    with open("test-data/config_1.yaml") as f:
        config = yaml.safe_load(f)
    # first input is the model_building input containing Mpro-IBM0078
    config["inputs"][0][Constants.CONFIG_EXCLUDE] = list(exclude)
    with open(version_dir / "config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)

    return tmp_path, version_dir


def _run_collator(working_dir: Path, force_superseded: bool):
    c = Collator(str(working_dir), no_validate_configs=True, force_superseded=force_superseded)
    meta = c.validate()
    assert meta is not None and len(c.logger.errors) == 0, c.logger.errors
    c.run(meta)


def test_force_superseded_status_and_files(tmp_path):
    working_dir, version_dir = _setup_working_dir(tmp_path, exclude=[EXCLUDED_DATASET])

    # --- upload_1: excluded dataset absent, everything else is "new" ---
    _run_collator(working_dir, force_superseded=False)

    meta_1 = _read_meta(version_dir / "upload_1")
    xtals_1 = meta_1[Constants.META_XTALS]
    assert EXCLUDED_DATASET not in xtals_1
    for dataset in UNCHANGED_DATASETS:
        assert xtals_1[dataset][Constants.META_STATUS] == Constants.META_STATUS_NEW

    # --- upload_2: drop the exclusion and force superseded ---
    # remove the exclude from the config, then create the new (empty) upload_2 dir
    with open(version_dir / "config.yaml") as f:
        config = yaml.safe_load(f)
    del config["inputs"][0][Constants.CONFIG_EXCLUDE]
    with open(version_dir / "config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)
    (version_dir / "upload_2").mkdir()

    _run_collator(working_dir, force_superseded=True)

    meta_2 = _read_meta(version_dir / "upload_2")
    xtals_2 = meta_2[Constants.META_XTALS]

    # the re-introduced dataset is brand new
    assert xtals_2[EXCLUDED_DATASET][Constants.META_STATUS] == Constants.META_STATUS_NEW

    # datasets whose PDB did not change are forced to superseded ...
    for dataset in UNCHANGED_DATASETS:
        assert xtals_2[dataset][Constants.META_STATUS] == Constants.META_STATUS_SUPERSEDED

        # ... and their files are re-copied into upload_2 with metadata pointing there
        pdb_entry = xtals_2[dataset][Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        pdb_file = pdb_entry[Constants.META_FILE]
        assert pdb_file.startswith("upload_2/"), pdb_file
        assert (version_dir / pdb_file).is_file()

    # the new dataset's files also live under upload_2
    new_pdb_file = xtals_2[EXCLUDED_DATASET][Constants.META_XTAL_FILES][Constants.META_XTAL_PDB][Constants.META_FILE]
    assert new_pdb_file.startswith("upload_2/"), new_pdb_file
    assert (version_dir / new_pdb_file).is_file()
