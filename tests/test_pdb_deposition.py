import textwrap
from pathlib import Path

import pytest
from gemmi import cif

from pdbdepo.pdb_deposition import (
    merge_mmcifgen_into_structure,
    read_cmpd_codes,
    read_fragalysis_csv,
    substitute_tokens,
)


# ---------------------------------------------------------------------------
# substitute_tokens — pure unit tests
# ---------------------------------------------------------------------------


def test_substitute_tokens_compound_and_crystal():
    result = substitute_tokens('$CompoundCode in $CrystalName', 'XTAL001', 'ABC123', None, '')
    assert result == 'ABC123 in XTAL001'


def test_substitute_tokens_pose_id_single():
    result = substitute_tokens('keywords, $PoseID', 'XTAL001', 'ABC', None, 'XTAL001a')
    assert result == 'keywords, XTAL001a'


def test_substitute_tokens_pose_id_multiple():
    result = substitute_tokens('keywords, $PoseID', 'XTAL001', 'ABC', None, 'XTAL001a, XTAL001b')
    assert result == 'keywords, XTAL001a, XTAL001b'


def test_substitute_tokens_external_code_correct_column():
    # ext_codes[0] -> $ExternalCode2, ext_codes[1] -> $ExternalCode3
    result = substitute_tokens('id $ExternalCode3', 'XTAL001', 'ABC', ['OB-001', 'OB-002'], '')
    assert result == 'id OB-002'


def test_substitute_tokens_erases_unsubstituted_external_codes():
    result = substitute_tokens('($ExternalCode3)', 'XTAL001', 'ABC', None, '')
    assert result == '()'


def test_substitute_tokens_no_ext_codes_erases_all():
    result = substitute_tokens('$ExternalCode2 $ExternalCode5', 'X', 'C', None, '')
    assert result == ' '


def test_substitute_tokens_no_tokens_unchanged():
    result = substitute_tokens('no tokens here', 'XTAL001', 'ABC', None, '')
    assert result == 'no tokens here'


def test_substitute_tokens_empty_pose_id_leaves_no_token():
    result = substitute_tokens('keys, $PoseID', 'XTAL001', 'ABC', None, '')
    assert '$PoseID' not in result
    assert result == 'keys, '


# ---------------------------------------------------------------------------
# read_fragalysis_csv
# ---------------------------------------------------------------------------


def test_read_fragalysis_csv_happy_path(tmp_path):
    csv = tmp_path / 'meta.csv'
    csv.write_text(
        textwrap.dedent(
            """\
        Code,Long code,Experiment code,Compound code
        A71EV2A-x0836a,A71EV2A-x0836_A_301_0_v1,A71EV2A-x0836,Z123
        A71EV2A-x0836b,A71EV2A-x0836_A_302_0_v1,A71EV2A-x0836,Z123
        A71EV2A-x0001a,A71EV2A-x0001_A_1_0_v1,A71EV2A-x0001,Z456
        """
        )
    )
    result = read_fragalysis_csv(str(csv))
    assert result == {
        'A71EV2A-x0836': ['A71EV2A-x0836a', 'A71EV2A-x0836b'],
        'A71EV2A-x0001': ['A71EV2A-x0001a'],
    }


def test_read_fragalysis_csv_quoted_fields(tmp_path):
    """Quoted CSV fields (as produced by the real Fragalysis download) must be parsed correctly."""
    csv = tmp_path / 'meta.csv'
    csv.write_text(
        textwrap.dedent(
            """\
        Code,Long code,Experiment code,Compound code
        "A71EV2A-x0836a","A71EV2A-x0836_A_301_0_v1","A71EV2A-x0836","Z123"
        "A71EV2A-x0836b","A71EV2A-x0836_A_302_0_v1","A71EV2A-x0836","Z123"
        """
        )
    )
    result = read_fragalysis_csv(str(csv))
    assert result == {'A71EV2A-x0836': ['A71EV2A-x0836a', 'A71EV2A-x0836b']}


def test_read_fragalysis_csv_none_filename():
    assert read_fragalysis_csv(None) == {}


def test_read_fragalysis_csv_skips_short_rows(tmp_path):
    csv = tmp_path / 'meta.csv'
    csv.write_text('Code,Long code,Experiment code\nA-x0001a,long\n')  # only 2 columns on data row
    result = read_fragalysis_csv(str(csv))
    assert result == {}


# ---------------------------------------------------------------------------
# read_cmpd_codes
# ---------------------------------------------------------------------------


def test_read_cmpd_codes_happy_path(tmp_path):
    csv = tmp_path / 'codes.csv'
    csv.write_text(
        textwrap.dedent(
            """\
        CrystalName,CompoundCode,OpenBindId
        XTAL001,Z001,OB-001
        XTAL002,Z002,OB-002
        """
        )
    )
    result = read_cmpd_codes(str(csv))
    assert result == {'XTAL001': ['Z001', 'OB-001'], 'XTAL002': ['Z002', 'OB-002']}


def test_read_cmpd_codes_none_filename():
    assert read_cmpd_codes(None) == {}


# ---------------------------------------------------------------------------
# merge_mmcifgen_into_structure — CIF-level tests (in-memory gemmi objects)
# ---------------------------------------------------------------------------


def _make_mmcifgen_doc():
    """Return a minimal mmcif-gen CIF document with a _struct.title loop containing tokens."""
    doc = cif.Document()
    block = doc.add_new_block('mmcifgen')
    loop = block.init_loop('', ['_struct.entry_id', '_struct.title'])
    loop.add_row(['INVID', 'structure of $CrystalName with $CompoundCode ($ExternalCode3)'])
    return doc


def _make_struct_block():
    """Return an empty destination CIF block."""
    doc = cif.Document()
    return doc.add_new_block('output')


def test_merge_mmcifgen_title_token_substitution():
    mmcifgen_block = _make_mmcifgen_doc()[0]
    struct_block = _make_struct_block()

    merge_mmcifgen_into_structure(
        struct_block,
        mmcifgen_block,
        xtal_name='XTAL042',
        cmpd_code='Z999',
        cmpd_codes_dict={'XTAL042': ['COL2', 'OB-042']},
        pose_ids_dict={},
    )

    title_item = next(item for item in struct_block if item.loop is not None and '_struct.title' in item.loop.tags)
    title = title_item.loop.values[1]
    assert 'XTAL042' in title
    assert 'Z999' in title
    assert 'OB-042' in title
    assert '$' not in title


def test_merge_mmcifgen_title_erases_unused_external_code():
    mmcifgen_block = _make_mmcifgen_doc()[0]
    struct_block = _make_struct_block()

    merge_mmcifgen_into_structure(
        struct_block,
        mmcifgen_block,
        xtal_name='XTAL042',
        cmpd_code='Z999',
        cmpd_codes_dict={},  # no external codes for this crystal
        pose_ids_dict={},
    )

    title_item = next(item for item in struct_block if item.loop is not None and '_struct.title' in item.loop.tags)
    title = title_item.loop.values[1]
    assert '$ExternalCode3' not in title
    assert '$' not in title


def test_merge_mmcifgen_keywords_loop_pose_id():
    """A keywords loop containing $PoseID in any value has it substituted."""
    doc = cif.Document()
    block = doc.add_new_block('mmcifgen')
    loop = block.init_loop(
        '', ['_struct_keywords.entry_id', '_struct_keywords.text', '_struct_keywords.pdbx_keywords']
    )
    loop.add_row(['INVID', 'Diamond Light Source, $PoseID', 'VIRAL PROTEIN'])

    struct_block = _make_struct_block()
    merge_mmcifgen_into_structure(
        struct_block,
        block,
        xtal_name='XTAL042',
        cmpd_code='Z999',
        cmpd_codes_dict={},
        pose_ids_dict={'XTAL042': ['XTAL042a', 'XTAL042b']},
    )

    kw_item = next(
        item for item in struct_block if item.loop is not None and '_struct_keywords.text' in item.loop.tags
    )
    idx = kw_item.loop.tags.index('_struct_keywords.text')
    text_value = kw_item.loop.values[idx]
    assert 'XTAL042a, XTAL042b' in text_value
    assert '$PoseID' not in text_value


def test_merge_mmcifgen_keywords_pair_pose_id():
    """A pair item whose value contains $PoseID has it substituted."""
    doc = cif.Document()
    block = doc.add_new_block('mmcifgen')
    block.set_pair('_struct_keywords.text', 'Diamond Light Source, $PoseID')

    struct_block = _make_struct_block()
    merge_mmcifgen_into_structure(
        struct_block,
        block,
        xtal_name='XTAL042',
        cmpd_code='Z999',
        cmpd_codes_dict={},
        pose_ids_dict={'XTAL042': ['XTAL042a']},
    )

    kw_pair = next(item for item in struct_block if item.pair is not None and '_struct_keywords.text' in item.pair[0])
    assert 'XTAL042a' in kw_pair.pair[1]
    assert '$PoseID' not in kw_pair.pair[1]


def test_merge_mmcifgen_passthrough_item_added_unchanged():
    """Items with no tokens are added to the destination block as-is."""
    doc = cif.Document()
    block = doc.add_new_block('mmcifgen')
    block.set_pair('_exptl.method', 'X-RAY DIFFRACTION')

    struct_block = _make_struct_block()
    merge_mmcifgen_into_structure(struct_block, block, 'XTAL', 'Z1', {}, {})

    exptl = next((item for item in struct_block if item.pair is not None and item.pair[0] == '_exptl.method'), None)
    assert exptl is not None
    assert exptl.pair[1] == 'X-RAY DIFFRACTION'


def test_merge_mmcifgen_no_pose_ids_for_crystal():
    """A crystal not present in pose_ids_dict produces an empty substitution — no crash, no token."""
    doc = cif.Document()
    block = doc.add_new_block('mmcifgen')
    loop = block.init_loop('', ['_struct_keywords.entry_id', '_struct_keywords.text'])
    loop.add_row(['INVID', 'keywords, $PoseID'])

    struct_block = _make_struct_block()
    merge_mmcifgen_into_structure(struct_block, block, 'MISSING_XTAL', 'Z1', {}, {})

    kw_item = next(
        item for item in struct_block if item.loop is not None and '_struct_keywords.text' in item.loop.tags
    )
    idx = kw_item.loop.tags.index('_struct_keywords.text')
    assert '$PoseID' not in kw_item.loop.values[idx]
