import textwrap
from pathlib import Path

import pandas as pd
import pytest
from gemmi import cif

from pdbdepo import pdb_deposition
from pdbdepo.pdb_deposition import (
    filter_excluded,
    merge_mmcifgen_into_structure,
    read_cmpd_codes,
    read_fragalysis_csv,
    rename_beamlines,
    substitute_tokens,
    validate_sequences,
)
from xchemalign.utils import Constants


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


# ---------------------------------------------------------------------------
# rename_beamlines — CIF-level tests (in-memory gemmi objects)
# ---------------------------------------------------------------------------


def _make_diffrn_source_block(*beamlines):
    """Return a block with a _diffrn_source loop, one row per beamline name given."""
    doc = cif.Document()
    block = doc.add_new_block('TEST')
    loop = block.init_loop(
        '',
        [
            '_diffrn_source.source',
            '_diffrn_source.type',
            '_diffrn_source.pdbx_synchrotron_site',
            '_diffrn_source.pdbx_synchrotron_beamline',
            '_diffrn_source.diffrn_id',
        ],
    )
    for i, beamline in enumerate(beamlines):
        loop.add_row(['SYNCHROTRON', cif.quote('DIAMOND BEAMLINE ' + beamline), 'DIAMOND', beamline, str(i + 1)])
    return doc, block


def _beamline_values(block):
    """Return (type, pdbx_synchrotron_beamline) unquoted values, one tuple per row."""
    types = block.find_values('_diffrn_source.type')
    beamlines = block.find_values('_diffrn_source.pdbx_synchrotron_beamline')
    return [(types.str(i), beamlines.str(i)) for i in range(len(types))]


def test_rename_beamlines_i02_2_to_vmxi():
    doc, block = _make_diffrn_source_block('I02-2')

    changes = rename_beamlines(block)

    assert _beamline_values(block) == [('DIAMOND BEAMLINE VMXi', 'VMXi')]
    assert len(changes) == 2


def test_rename_beamlines_i02_1_to_vmxm():
    doc, block = _make_diffrn_source_block('I02-1')

    rename_beamlines(block)

    assert _beamline_values(block) == [('DIAMOND BEAMLINE VMXm', 'VMXm')]


@pytest.mark.parametrize('beamline', ['I02', 'I03', 'I04', 'I04-1', 'I23', 'I24', 'VMXi', 'VMXm'])
def test_rename_beamlines_leaves_wwpdb_names_untouched(beamline):
    """Names already in the wwPDB enumeration must pass through unchanged — notably I02, which
    must not be caught by the I02-1/I02-2 rules."""
    doc, block = _make_diffrn_source_block(beamline)

    changes = rename_beamlines(block)

    assert changes == []
    assert _beamline_values(block) == [('DIAMOND BEAMLINE ' + beamline, beamline)]


def test_rename_beamlines_no_diffrn_source():
    """A block with no _diffrn_source at all is a no-op, not an error."""
    block = _make_struct_block()

    assert rename_beamlines(block) == []


def test_rename_beamlines_pair_form():
    """_diffrn_source written as pairs rather than a loop is renamed too."""
    doc = cif.Document()
    block = doc.add_new_block('TEST')
    block.set_pair('_diffrn_source.type', cif.quote('DIAMOND BEAMLINE I02-2'))
    block.set_pair('_diffrn_source.pdbx_synchrotron_beamline', 'I02-2')

    rename_beamlines(block)

    assert _beamline_values(block) == [('DIAMOND BEAMLINE VMXi', 'VMXi')]


def test_rename_beamlines_multi_row_only_renames_matching():
    doc, block = _make_diffrn_source_block('I04-1', 'I02-2')

    rename_beamlines(block)

    assert _beamline_values(block) == [
        ('DIAMOND BEAMLINE I04-1', 'I04-1'),
        ('DIAMOND BEAMLINE VMXi', 'VMXi'),
    ]


def test_rename_beamlines_preserves_quoting():
    """The rewritten values must be requoted: the type string has spaces, the bare name does not."""
    doc, block = _make_diffrn_source_block('I02-2')

    rename_beamlines(block)

    written = doc.as_string()
    assert "'DIAMOND BEAMLINE VMXi' DIAMOND VMXi" in written
    assert 'I02-2' not in written


# ---------------------------------------------------------------------------
# filter_excluded
# ---------------------------------------------------------------------------


def _soakdb_df(*xtal_names):
    return pd.DataFrame(
        {Constants.SOAKDB_XTAL_NAME: list(xtal_names), 'RefinementOutcome': ['5 - Deposition'] * len(xtal_names)}
    )


def test_filter_excluded_drops_listed_crystals():
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002', 'XTAL-x0003')
    result = filter_excluded(df, {'exclude': ['XTAL-x0002']})
    assert list(result[Constants.SOAKDB_XTAL_NAME]) == ['XTAL-x0001', 'XTAL-x0003']


def test_filter_excluded_drops_several_crystals():
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002', 'XTAL-x0003')
    result = filter_excluded(df, {'exclude': ['XTAL-x0001', 'XTAL-x0003']})
    assert list(result[Constants.SOAKDB_XTAL_NAME]) == ['XTAL-x0002']


def test_filter_excluded_no_exclude_section():
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002')
    result = filter_excluded(df, {'dir': 'somewhere'})
    assert list(result[Constants.SOAKDB_XTAL_NAME]) == ['XTAL-x0001', 'XTAL-x0002']


def test_filter_excluded_empty_exclude_tag():
    # 'exclude:' present in the YAML but with no values parses as None
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002')
    result = filter_excluded(df, {'exclude': None})
    assert list(result[Constants.SOAKDB_XTAL_NAME]) == ['XTAL-x0001', 'XTAL-x0002']


def test_filter_excluded_name_not_in_dataframe():
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002')
    result = filter_excluded(df, {'exclude': ['XTAL-x9999']})
    assert list(result[Constants.SOAKDB_XTAL_NAME]) == ['XTAL-x0001', 'XTAL-x0002']


def test_filter_excluded_leaves_other_columns_intact():
    df = _soakdb_df('XTAL-x0001', 'XTAL-x0002')
    result = filter_excluded(df, {'exclude': ['XTAL-x0001']})
    assert list(result.columns) == list(df.columns)
    assert result.iloc[0]['RefinementOutcome'] == '5 - Deposition'


def test_excluded_crystal_is_not_sequence_checked(monkeypatch, tmp_path):
    # the reported bug: a crystal in the exclude section still failed the issue #94 sequence check,
    # which is fatal, so the user had no way of getting past it
    monkeypatch.setattr(pdb_deposition.sequence_check, 'collect_candidates', lambda *args: {})
    monkeypatch.setattr(pdb_deposition.sequence_check, 'read_structure', lambda **kwargs: 'structure')

    checked = []

    def fake_check_sequences(struc, seq_dict):
        checked.append(seq_dict)
        return ['chain A: declared D but model has LEU 33']

    monkeypatch.setattr(pdb_deposition.sequence_check, 'check_sequences', fake_check_sequences)
    monkeypatch.setattr(pdb_deposition.sequence_check, 'suggest_sequence', lambda struc, candidates: None)

    df = pd.DataFrame(
        {
            Constants.SOAKDB_XTAL_NAME: ['XTAL-x0001'],
            Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST: ['None'],
            Constants.SOAKDB_COL_PDB: [str(tmp_path / 'XTAL-x0001.pdb')],
        }
    )
    input_config = {'exclude': ['XTAL-x0001']}
    default_seq = {'A': ('A', 'MNPQRSTVW')}

    # without the filter the failing crystal is checked, and is fatal
    with pytest.raises(SystemExit):
        validate_sequences(tmp_path, df, input_config, default_seq, {})
    assert checked == [default_seq]

    # once it is excluded it is not checked at all
    checked.clear()
    filtered = filter_excluded(df, input_config)
    validate_sequences(tmp_path, filtered, input_config, default_seq, {})
    assert not checked
