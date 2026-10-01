# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0

"""Unit tests for Aligner._build_aligned_files, in particular the matching of event map
metadata to aligned ligands when the altlocs change between incremental alignments."""

from types import SimpleNamespace
from unittest.mock import MagicMock

from xchemalign import utils
from xchemalign.aligner import Aligner

C = utils.Constants

SITE = "1"


def _aligner():
    a = Aligner.__new__(Aligner)
    a.logger = MagicMock()
    return a


def _version_output(tag):
    """The LNA output for one ligand observation, with a single site."""
    return SimpleNamespace(
        aligned_structures={SITE: f"{tag}.pdb"},
        aligned_artefacts={SITE: f"{tag}_artefacts.pdb"},
        aligned_event_maps={SITE: f"{tag}_event.ccp4"},
        aligned_xmaps={SITE: f"{tag}_xmap.ccp4"},
        aligned_diff_maps={SITE: f"{tag}_diff.ccp4"},
        aligned_event_maps_crystallographic={SITE: f"{tag}_event_x.ccp4"},
        aligned_xmaps_crystallographic={SITE: f"{tag}_xmap_x.ccp4"},
        aligned_diff_maps_crystallographic={SITE: f"{tag}_diff_x.ccp4"},
    )


def _event(chain, res, altloc, with_file=True):
    e = {C.META_PROT_CHAIN: chain, C.META_PROT_RES: res, C.META_PROT_ALTLOC: altloc}
    if with_file:
        e[C.META_FILE] = f"{chain}_{res}_{altloc}.ccp4"
    return e


def _has_event_map(aligned_version):
    return C.META_AIGNED_EVENT_MAP in aligned_version[SITE]


def test_altloc_added_between_versions():
    """A ligand previously modelled as one altloc is now modelled as two. Only the original altloc
    has event map metadata. This used to exit with 'Unexpected number of event maps'."""
    dataset_output = {"A": {"201": {"A": {"v1": _version_output("a")}, "B": {"v1": _version_output("b")}}}}
    events = [_event("A", 201, "A")]

    out = _aligner()._build_aligned_files("xtal-1", dataset_output, events)

    assert _has_event_map(out["A"]["201"]["A"]["v1"])
    assert not _has_event_map(out["A"]["201"]["B"]["v1"])


def test_event_map_matched_by_key_not_position():
    """The events list is in a different order to the alignment output, and has an extra entry for a
    ligand with no alignment output. Each altloc must still get its own event map status."""
    dataset_output = {"A": {"201": {"A": {"v1": _version_output("a")}, "B": {"v1": _version_output("b")}}}}
    events = [
        _event("A", 999, "A"),
        _event("A", 201, "B", with_file=False),
        _event("A", 201, "A"),
    ]

    out = _aligner()._build_aligned_files("xtal-1", dataset_output, events)

    assert _has_event_map(out["A"]["201"]["A"]["v1"])
    assert not _has_event_map(out["A"]["201"]["B"]["v1"])


def test_null_altloc_matches_string_zero():
    """LNA uses '\\0' for no altloc, whereas the collator metadata records it as '0'."""
    dataset_output = {"A": {"201": {"\0": {"v1": _version_output("a")}}}}
    events = [_event("A", 201, "0")]

    out = _aligner()._build_aligned_files("xtal-1", dataset_output, events)

    assert _has_event_map(out["A"]["201"]["\0"]["v1"])


def test_event_map_not_matched_across_chains():
    """The same residue and altloc in different chains must not share event map metadata."""
    dataset_output = {
        "A": {"201": {"A": {"v1": _version_output("a")}}},
        "B": {"201": {"A": {"v1": _version_output("b")}}},
    }
    events = [_event("B", 201, "A")]

    out = _aligner()._build_aligned_files("xtal-1", dataset_output, events)

    assert not _has_event_map(out["A"]["201"]["A"]["v1"])
    assert _has_event_map(out["B"]["201"]["A"]["v1"])


def test_missing_event_map_metadata_warns_and_continues():
    dataset_output = {"A": {"201": {"A": {"v1": _version_output("a")}}}}
    aligner = _aligner()

    out = aligner._build_aligned_files("xtal-1", dataset_output, [])

    assert not _has_event_map(out["A"]["201"]["A"]["v1"])
    assert out["A"]["201"]["A"]["v1"][SITE][C.META_AIGNED_STRUCTURE] == "a.pdb"
    aligner.logger.warn.assert_called_once()
    assert "xtal-1" in aligner.logger.warn.call_args[0][0]
