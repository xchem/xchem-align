# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0

"""Unit tests for the compound-curation additions in xchemalign.uploader:

* extracting incoming compounds from meta_aligner.yaml for the validation call;
* sending them (and any completed curation file) as JSON;
* reacting to a ``compound_conflicts`` validation response by writing the
  server-generated spreadsheet and halting the upload.
"""

import base64
from unittest.mock import MagicMock

import pytest

from xchemalign import utils
from xchemalign.uploader import (
    CurationRequiredError,
    TarballSource,
    XCADataUpload,
)

C = utils.Constants


def _meta(*ligands_per_crystal):
    """Build a minimal meta_aligner.yaml dict from lists of ligand dicts."""
    crystals = {}
    for i, ligands in enumerate(ligands_per_crystal):
        crystals[f"Xtal-{i}"] = {
            C.META_XTAL_FILES: {
                C.META_XTAL_CIF: {
                    C.META_LIGANDS: {f"LIG{j}": lig for j, lig in enumerate(ligands)},
                },
            },
        }
    return {C.META_XTALS: crystals}


def _source(meta):
    src = TarballSource.__new__(TarballSource)
    src._meta = meta
    src._config = {"target_name": "MyTarget"}
    src._target_name = "MyTarget"
    src._data_version = "1.0"
    src._upload_version = "1"
    return src


def test_extract_compounds_basic():
    src = _source(
        _meta(
            [{C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1", C.META_MODELED_SMILES_CANON: "CCO"}],
        )
    )
    compounds = src._extract_compounds()
    assert compounds == [{"smiles": "CCO", "compound_code": "Z1", "modeled_smiles_canon": "CCO"}]


def test_extract_compounds_dedupes_identical_across_crystals():
    lig = {C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1"}
    src = _source(_meta([lig], [dict(lig)]))
    assert src._extract_compounds() == [{"smiles": "CCO", "compound_code": "Z1"}]


def test_extract_compounds_keeps_distinct_same_smiles_different_code():
    src = _source(
        _meta(
            [
                {C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1"},
                {C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z2"},
            ],
        )
    )
    assert len(src._extract_compounds()) == 2


def test_extract_compounds_skips_ligand_without_smiles():
    src = _source(_meta([{C.META_CMPD_CODE: "Z1"}]))
    assert not src._extract_compounds()


def test_extract_compounds_omits_missing_and_null_fields():
    src = _source(_meta([{C.META_SMILES: "CCO", C.META_CMPD_CODE: None}]))
    assert src._extract_compounds() == [{"smiles": "CCO"}]


def test_validation_payload_includes_compounds():
    src = _source(_meta([{C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1"}]))
    payload = src.get_validation_payload()
    assert payload["compounds"] == [{"smiles": "CCO", "compound_code": "Z1"}]
    assert payload["target_name"] == "MyTarget"


def _uploader(data_source, curation_file=None):
    up = XCADataUpload(
        url="https://example.com",
        data_source=data_source,
        proposal="lb-1",
        curation_file=str(curation_file) if curation_file else None,
    )
    up._session = MagicMock()
    return up


def _resp(json_data, ok=True, url="https://example.com/api/validate"):
    r = MagicMock()
    r.ok = ok
    r.url = url
    r.json.return_value = json_data
    return r


def test_validate_sends_json_payload():
    src = _source(_meta([{C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1"}]))
    up = _uploader(src)
    up._session.post.return_value = _resp({"success": True, "message": []})

    up._validate()

    _, kwargs = up._session.post.call_args
    assert "json" in kwargs and "data" not in kwargs
    assert kwargs["json"]["compounds"] == [{"smiles": "CCO", "compound_code": "Z1"}]
    assert kwargs["json"]["target_access_string"] == "lb-1"


def test_validate_conflicts_writes_xlsx_and_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    src = _source(_meta([{C.META_SMILES: "CCO", C.META_CMPD_CODE: "Z1"}]))
    up = _uploader(src)
    xlsx_bytes = b"PK\x03\x04 fake xlsx"
    up._session.post.return_value = _resp(
        {
            "success": False,
            "message": ["1 compound(s) need review before upload"],
            "compound_conflicts": [{"inchi_key": "AAA", "status": "conflict"}],
            "curation_file": base64.b64encode(xlsx_bytes).decode("ascii"),
            "curation_filename": "MyTarget_curation.xlsx",
        }
    )

    with pytest.raises(CurationRequiredError) as exc:
        up._validate()

    written = tmp_path / "MyTarget_curation.xlsx"
    assert written.read_bytes() == xlsx_bytes
    assert exc.value.curation_path == written


def test_validate_with_curation_file_sends_it_encoded(tmp_path):
    curation = tmp_path / "done.xlsx"
    curation.write_bytes(b"PK\x03\x04 resolved")
    src = _source(_meta([{C.META_SMILES: "CCO"}]))
    up = _uploader(src, curation_file=curation)
    up._session.post.return_value = _resp({"success": True, "message": []})

    up._validate()

    _, kwargs = up._session.post.call_args
    sent = base64.b64decode(kwargs["json"]["curation_file"])
    assert sent == b"PK\x03\x04 resolved"


def test_validate_no_conflicts_does_not_write(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    src = _source(_meta([{C.META_SMILES: "CCO"}]))
    up = _uploader(src)
    up._session.post.return_value = _resp({"success": True, "message": []})

    up._validate()

    assert not list(tmp_path.iterdir())
