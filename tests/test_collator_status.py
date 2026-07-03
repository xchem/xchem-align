from xchemalign.collator import Collator
from xchemalign.utils import Constants


def _make_collator(force_superseded=False, rejected_xtals=None):
    # bypass the heavy __init__ (config/working dir); exercise only _get_xtal_status
    c = Collator.__new__(Collator)
    c.force_superseded = force_superseded
    c.rejected_xtals = rejected_xtals or set()
    return c


def _xtal_data(sha256):
    return {Constants.META_XTAL_FILES: {Constants.META_XTAL_PDB: {Constants.META_SHA256: sha256}}}


def test_status_unchanged_when_pdb_identical():
    c = _make_collator()
    old = _xtal_data("aaa")
    new = _xtal_data("aaa")
    assert c._get_xtal_status("xtal-1", old, new) == Constants.META_STATUS_UNCHANGED


def test_status_superseded_when_pdb_differs():
    c = _make_collator()
    old = _xtal_data("aaa")
    new = _xtal_data("bbb")
    assert c._get_xtal_status("xtal-1", old, new) == Constants.META_STATUS_SUPERSEDED


def test_force_superseded_promotes_unchanged():
    c = _make_collator(force_superseded=True)
    old = _xtal_data("aaa")
    new = _xtal_data("aaa")
    assert c._get_xtal_status("xtal-1", old, new) == Constants.META_STATUS_SUPERSEDED


def test_rejected_xtal_stays_deprecated_regardless_of_force():
    c = _make_collator(force_superseded=True, rejected_xtals={"xtal-1"})
    old = _xtal_data("aaa")
    new = _xtal_data("aaa")
    assert c._get_xtal_status("xtal-1", old, new) == Constants.META_STATUS_DEPRECATED
