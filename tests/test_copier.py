from pathlib import Path

import pytest

from xchemalign import utils
from xchemalign.copier import Copier


def make_copier(tmp_path, sequence_file_paths):
    """A Copier wired up just enough to exercise copy_sequence_files, which touches no database."""
    return Copier(
        base_path=tmp_path / 'base',
        input_path=Path('visit/dir'),
        output_path=tmp_path / 'out',
        soakdb_file_path=Path(utils.Constants.DEFAULT_SOAKDB_PATH),
        panddas_file_paths=[],
        ref_datasets=[],
        statuses=None,
        sequence_file_paths=sequence_file_paths,
    )


def write_sequences(tmp_path, *names, seq_dir='processing/analysis/sequences'):
    d = tmp_path / 'base' / 'visit/dir' / seq_dir
    d.mkdir(parents=True, exist_ok=True)
    for name in names:
        (d / name).write_text('> A A\nACDEFGHIK\n')
    return [Path(seq_dir) / name for name in names]


def test_copies_the_default_and_the_variants(tmp_path):
    paths = write_sequences(tmp_path, 'default.fa', '2-chain.fa')
    c = make_copier(tmp_path, paths)

    assert c.copy_sequence_files() == 2
    out = tmp_path / 'out/visit/dir/processing/analysis/sequences'
    assert (out / 'default.fa').is_file()
    assert (out / '2-chain.fa').is_file()
    assert (out / 'default.fa').read_text() == '> A A\nACDEFGHIK\n'
    assert not c.errors


def test_missing_file_warns_and_the_others_are_still_copied(tmp_path):
    paths = write_sequences(tmp_path, 'default.fa')
    paths.append(Path('processing/analysis/sequences/gone.fa'))
    c = make_copier(tmp_path, paths)

    assert c.copy_sequence_files() == 1
    assert (tmp_path / 'out/visit/dir/processing/analysis/sequences/default.fa').is_file()
    assert len(c.warnings) == 1
    assert 'gone.fa' in c.warnings[0]
    # a missing sequence file must not stop the copy, only pdb_deposition treats this as fatal
    assert not c.errors


def test_no_sequences_declared_copies_nothing(tmp_path):
    c = make_copier(tmp_path, [])
    assert c.copy_sequence_files() == 0
    assert not c.warnings
    assert not c.errors


def test_sequences_are_copied_to_the_path_read_sequences_expects(tmp_path):
    """The copier and read_sequences must agree on where the files live, or deposition breaks."""
    input_yaml = {
        utils.Constants.CONFIG_DIR: 'visit/dir',
        utils.Constants.CONFIG_TYPE: utils.Constants.CONFIG_TYPE_MODEL_BUILDING,
        utils.Constants.CONFIG_SEQUENCES: {
            utils.Constants.CONFIG_DIR: 'processing/sequences',
            utils.Constants.CONFIG_DEFAULT: 'default.fa',
            utils.Constants.CONFIG_VARIANTS: [
                {utils.Constants.CONFIG_SEQUENCE: '2-chain.fa', utils.Constants.CONFIG_CRYSTALS: ['xtal-1']}
            ],
        },
    }
    write_sequences(tmp_path, 'default.fa', '2-chain.fa', seq_dir='processing/sequences')

    c = make_copier(tmp_path, utils.sequence_file_paths(input_yaml))
    assert c.copy_sequence_files() == 2

    # now read them back out of the copied tree the way pdb_deposition would
    default_seq, variants = utils.read_sequences(tmp_path / 'out', input_yaml)
    assert default_seq == {'A': ('A', 'ACDEFGHIK')}
    assert sorted(variants) == ['xtal-1']
