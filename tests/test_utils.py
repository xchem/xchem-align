from xchemalign import utils
from xchemalign.utils import Constants


def write_fasta(tmp_path, name, contents):
    p = tmp_path / name
    p.write_text(contents)
    return p


def test_read_fasta_single_entity(tmp_path):
    p = write_fasta(tmp_path, 'default.fa', '> A A\nACDEFGHIK\n')
    assert utils.read_fasta(p) == {'A': ('A', 'ACDEFGHIK')}


def test_read_fasta_one_entity_two_chains(tmp_path):
    p = write_fasta(tmp_path, '2-chain.fa', '> A AB\nACDEF\nGHIK\n')
    assert utils.read_fasta(p) == {'A': ('A', 'ACDEFGHIK'), 'B': ('A', 'ACDEFGHIK')}


def test_read_fasta_keeps_the_sequence_of_every_entity(tmp_path):
    # the heterodimer example from the user guide: each record must keep its own sequence
    p = write_fasta(tmp_path, 'hetero.fa', '> A A\nACDEFGHIK\n> B B\nMNPQRSTVW\n')
    assert utils.read_fasta(p) == {'A': ('A', 'ACDEFGHIK'), 'B': ('B', 'MNPQRSTVW')}


def test_read_fasta_missing_file_is_not_fatal_when_asked(tmp_path):
    assert utils.read_fasta(tmp_path / 'nope.fa', fatal=False) is None


def test_read_sequences_without_a_sequences_section(tmp_path):
    input_yaml = {Constants.CONFIG_DIR: 'some/dir', Constants.CONFIG_TYPE: Constants.CONFIG_TYPE_MANUAL}
    assert utils.read_sequences(tmp_path, input_yaml) == (None, {})


def test_read_sequences_with_a_missing_fasta_is_not_fatal_when_asked(tmp_path):
    input_yaml = {
        Constants.CONFIG_DIR: 'some/dir',
        Constants.CONFIG_TYPE: Constants.CONFIG_TYPE_MODEL_BUILDING,
        Constants.CONFIG_SEQUENCES: {Constants.CONFIG_DEFAULT: 'default.fa'},
    }
    assert utils.read_sequences(tmp_path, input_yaml, fatal=False) == (None, {})


def test_read_sequences_reads_the_default_and_the_variants(tmp_path):
    seq_dir = tmp_path / 'some/dir/sequences'
    seq_dir.mkdir(parents=True)
    write_fasta(seq_dir, 'default.fa', '> A A\nACDEFGHIK\n')
    write_fasta(seq_dir, '2-chain.fa', '> A AB\nACDEFGHIK\n')
    input_yaml = {
        Constants.CONFIG_DIR: 'some/dir',
        Constants.CONFIG_TYPE: Constants.CONFIG_TYPE_MODEL_BUILDING,
        Constants.CONFIG_SEQUENCES: {
            Constants.CONFIG_DIR: 'sequences',
            Constants.CONFIG_DEFAULT: 'default.fa',
            Constants.CONFIG_VARIANTS: [
                {Constants.CONFIG_SEQUENCE: '2-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-1', 'xtal-2']}
            ],
        },
    }
    default_seq, variants = utils.read_sequences(tmp_path, input_yaml)
    assert default_seq == {'A': ('A', 'ACDEFGHIK')}
    assert sorted(variants) == ['xtal-1', 'xtal-2']
    assert sorted(variants['xtal-1']) == ['A', 'B']


def test_parse_compound_smiles():
    inputs = {'string1': [1, 1], 'string1;string2': [2, 1, 1], 'string1;string21 string22;string3': [3, 1, 2, 1]}
    for s, r in inputs.items():
        result = utils.parse_compound_smiles(s)
        assert len(result) == r[0]
        for i, v in enumerate(result):
            assert len(v) == r[i + 1]
    print('OK')
