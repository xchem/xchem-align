from pathlib import Path

import pytest
from gemmi import cif

from pdbdepo.pdb_deposition import validate_structure_cif_doc

CIF_DIR = Path(__file__).parent.parent / "test-data" / "cif"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _minimal_valid_doc():
    """Build a CIF document that passes all validation checks."""
    doc = cif.Document()
    block = doc.add_new_block('TEST')

    block.set_pair('_cell.length_a', '50.0')
    block.set_pair('_symmetry.space_group_name_H-M', 'P 21 21 21')

    loop = block.init_loop('', ['_atom_site.id', '_atom_site.label_entity_id'])
    loop.add_row(['1', '1'])

    loop = block.init_loop('', ['_entity.id', '_entity.type'])
    loop.add_row(['1', 'polymer'])

    loop = block.init_loop(
        '',
        [
            '_refine.entry_id',
            '_refine.ls_R_factor_R_work',
            '_refine.ls_R_factor_R_free',
            '_refine.ls_d_res_high',
            '_refine.ls_d_res_low',
            '_refine.B_iso_mean',
        ],
    )
    loop.add_row(['TEST', '0.20', '0.25', '1.5', '50.0', '20.0'])

    loop = block.init_loop('', ['_software.name', '_software.pdbx_ordinal'])
    loop.add_row(['REFMAC', '1'])

    return doc


# ---------------------------------------------------------------------------
# File-based tests (existing)
# ---------------------------------------------------------------------------


def test_disconnected_loop_raises_on_parse():
    """gemmi rejects a CIF where the same loop category appears in two separate, non-contiguous loops."""
    with pytest.raises(RuntimeError, match="duplicate tag"):
        cif.read(str(CIF_DIR / "disconnected_loop.cif"))


def test_duplicate_pair_raises_on_parse():
    """gemmi rejects a CIF with a duplicate pair tag."""
    with pytest.raises(RuntimeError, match="duplicate tag"):
        cif.read(str(CIF_DIR / "duplicate_pair.cif"))


def test_disconnected_pair_no_validation_issues():
    """A CIF with a disconnected pair (same category appearing both before and after a loop) parses without error and passes validation."""
    doc = cif.read(str(CIF_DIR / "disconnected_pair.cif"))
    issues = validate_structure_cif_doc(doc)
    assert issues == []


# ---------------------------------------------------------------------------
# In-memory validation tests
# ---------------------------------------------------------------------------


def test_minimal_valid_doc_passes():
    assert validate_structure_cif_doc(_minimal_valid_doc()) == []


def test_missing_required_loop_category_reported():
    doc = _minimal_valid_doc()
    # Remove the _refine loop
    block = doc[0]
    for item in block:
        if item.loop is not None and any(t.startswith('_refine.') for t in item.loop.tags):
            item.erase()
            break
    issues = validate_structure_cif_doc(doc)
    assert any('_refine' in i for i in issues)


def test_missing_required_pair_category_reported():
    doc = _minimal_valid_doc()
    block = doc[0]
    for item in block:
        if item.pair is not None and item.pair[0].startswith('_cell.'):
            item.erase()
            break
    issues = validate_structure_cif_doc(doc)
    assert any('_cell' in i for i in issues)


def _doc_with_rfactors(r_work, r_free):
    """Build a minimal valid doc with specific R-factor values."""
    doc = cif.Document()
    block = doc.add_new_block('TEST')
    block.set_pair('_cell.length_a', '50.0')
    block.set_pair('_symmetry.space_group_name_H-M', 'P 21 21 21')
    loop = block.init_loop('', ['_atom_site.id', '_atom_site.label_entity_id'])
    loop.add_row(['1', '1'])
    loop = block.init_loop('', ['_entity.id', '_entity.type'])
    loop.add_row(['1', 'polymer'])
    loop = block.init_loop(
        '',
        [
            '_refine.entry_id',
            '_refine.ls_R_factor_R_work',
            '_refine.ls_R_factor_R_free',
            '_refine.ls_d_res_high',
            '_refine.ls_d_res_low',
            '_refine.B_iso_mean',
        ],
    )
    loop.add_row(['TEST', str(r_work), str(r_free), '1.5', '50.0', '20.0'])
    loop = block.init_loop('', ['_software.name', '_software.pdbx_ordinal'])
    loop.add_row(['REFMAC', '1'])
    return doc


def test_rfactor_out_of_range_reported():
    doc = _doc_with_rfactors(r_work=1.5, r_free=0.25)  # R_work > 1.0 — invalid
    issues = validate_structure_cif_doc(doc)
    assert any('R_work' in i for i in issues)


def test_rfree_less_than_rwork_reported():
    doc = _doc_with_rfactors(r_work=0.30, r_free=0.20)  # R_free < R_work — invalid
    issues = validate_structure_cif_doc(doc)
    assert any('R_free' in i and 'R_work' in i for i in issues)


def test_unsubstituted_placeholder_in_title_reported():
    doc = _minimal_valid_doc()
    block = doc[0]
    loop = block.init_loop('', ['_struct.entry_id', '_struct.title'])
    loop.add_row(['TEST', 'structure of $CompoundCode ($CrystalName)'])
    issues = validate_structure_cif_doc(doc)
    assert any('placeholder' in i for i in issues)


def test_software_ordinal_sequential_passes():
    doc = _minimal_valid_doc()
    block = doc[0]
    # Replace single-entry software loop with a proper multi-entry sequential one
    for item in block:
        if item.loop is not None and '_software.pdbx_ordinal' in item.loop.tags:
            item.erase()
            break
    loop = block.init_loop('', ['_software.name', '_software.pdbx_ordinal'])
    loop.add_row(['REFMAC', '1'])
    loop.add_row(['CCP4', '2'])
    assert validate_structure_cif_doc(doc) == []


def test_software_ordinal_non_sequential_reported():
    doc = _minimal_valid_doc()
    block = doc[0]
    for item in block:
        if item.loop is not None and '_software.pdbx_ordinal' in item.loop.tags:
            item.erase()
            break
    loop = block.init_loop('', ['_software.name', '_software.pdbx_ordinal'])
    loop.add_row(['REFMAC', '1'])
    loop.add_row(['CCP4', '3'])  # gap — should be 2
    issues = validate_structure_cif_doc(doc)
    assert any('ordinal' in i for i in issues)


# ---------------------------------------------------------------------------
# Diamond beamline names (check 7)
# ---------------------------------------------------------------------------


def _doc_with_beamline(beamline):
    """Build a minimal valid doc carrying a _diffrn_source loop for the given beamline."""
    doc = _minimal_valid_doc()
    loop = doc[0].init_loop(
        '',
        ['_diffrn_source.type', '_diffrn_source.pdbx_synchrotron_beamline', '_diffrn_source.diffrn_id'],
    )
    loop.add_row([cif.quote('DIAMOND BEAMLINE ' + beamline), beamline, '1'])
    return doc


def test_ispyb_beamline_name_reported():
    issues = validate_structure_cif_doc(_doc_with_beamline('I02-2'))
    assert len(issues) == 1
    assert '_diffrn_source.type' in issues[0]
    assert 'I02-2' in issues[0]


def test_wwpdb_beamline_name_accepted():
    assert validate_structure_cif_doc(_doc_with_beamline('VMXi')) == []
    assert validate_structure_cif_doc(_doc_with_beamline('I04-1')) == []


def test_missing_diffrn_source_not_reported():
    """_diffrn_source is optional as far as this check is concerned."""
    assert validate_structure_cif_doc(_minimal_valid_doc()) == []
