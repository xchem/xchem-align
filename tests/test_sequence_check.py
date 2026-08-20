import copy

import gemmi
import pytest

from xchemalign import sequence_check

# chain A is complete (residues 1-302), chain B has an internal gap (48-50 are not modelled)
TWO_CHAIN_PDB = 'test-data/8e1y.pdb'
ONE_CHAIN_PDB = 'test-data/refine_7.split.bound-state.pdb'


def read(path):
    return sequence_check.read_structure(pdb_file=path)


def observed(struc, chain_name):
    return sequence_check.observed_sequence(sequence_check.polymer_chains(struc)[chain_name])


def seq_dict(struc, *chain_names):
    """Declare each named chain's own observed sequence, i.e. what a correct FASTA file would say."""
    return {name: ('A', observed(struc, name)) for name in chain_names}


def test_matching_sequence_has_no_issues():
    struc = read(ONE_CHAIN_PDB)
    assert sequence_check.check_sequences(struc, seq_dict(struc, 'A')) == []


def test_unmodelled_termini_are_not_reported():
    # the declared sequence is the full construct, the model only covers part of it, which is normal
    struc = read(ONE_CHAIN_PDB)
    declared = 'MGSSHHHHHHSSGLVPRGSH' + observed(struc, 'A') + 'GGSKLE'
    assert sequence_check.check_sequences(struc, {'A': ('A', declared)}) == []


def test_unmodelled_internal_residues_are_not_reported():
    # chain B of 8e1y does not model residues 48-50, so the declared sequence has three more residues
    struc = read(TWO_CHAIN_PDB)
    obs = observed(struc, 'B')
    declared = seq_dict(struc, 'A')
    declared['B'] = ('A', obs[:47] + 'AAA' + obs[47:])
    assert sequence_check.check_sequences(struc, declared) == []


def test_point_mutation_is_reported_with_the_residue_number():
    struc = read(ONE_CHAIN_PDB)
    obs = observed(struc, 'A')
    # residue 10 is at index 9 as this chain is numbered from 1
    declared = obs[:9] + ('W' if obs[9] != 'W' else 'Y') + obs[10:]
    issues = sequence_check.check_sequences(struc, {'A': ('A', declared)})
    assert len(issues) == 1
    assert 'chain A declares' in issues[0]
    assert issues[0].endswith(' 10')


def test_truncated_declared_sequence_is_reported():
    struc = read(ONE_CHAIN_PDB)
    issues = sequence_check.check_sequences(struc, {'A': ('A', observed(struc, 'A')[:-20])})
    assert issues
    assert all('does not account for' in issue for issue in issues)


def test_completely_wrong_sequence_is_reported():
    struc = read(ONE_CHAIN_PDB)
    declared = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR'
    assert sequence_check.check_sequences(struc, {'A': ('A', declared)})


def test_chain_in_model_but_not_declared():
    struc = read(TWO_CHAIN_PDB)
    issues = sequence_check.check_sequences(struc, seq_dict(struc, 'A'))
    assert issues == ['chain B is in the model but has no declared sequence']


def test_chain_declared_but_not_in_model():
    struc = read(ONE_CHAIN_PDB)
    declared = seq_dict(struc, 'A')
    declared['C'] = ('C', 'ACDEFGHIK')
    issues = sequence_check.check_sequences(struc, declared)
    assert issues == ['chain C has a declared sequence but is not in the model']


def test_both_chains_declared_has_no_issues():
    struc = read(TWO_CHAIN_PDB)
    assert sequence_check.check_sequences(struc, seq_dict(struc, 'A', 'B')) == []


@pytest.mark.parametrize('declared', [None, ''])
def test_empty_declared_sequence_is_reported(declared):
    struc = read(ONE_CHAIN_PDB)
    issues = sequence_check.check_sequences(struc, {'A': ('A', declared)})
    assert issues == ['chain A has an empty declared sequence']


def test_modified_residues_count_as_the_residue_they_derive_from():
    # a selenomethionine structure still matches a sequence that declares plain M
    struc = read(ONE_CHAIN_PDB)
    declared = observed(struc, 'A')
    modified = copy.deepcopy(struc)
    converted = 0
    for residue in sequence_check.polymer_chains(modified)['A']:
        if residue.name == 'MET':
            residue.name = 'MSE'
            converted += 1
    modified.setup_entities()
    assert converted > 0
    assert 'M' in declared
    assert sequence_check.check_sequences(modified, {'A': ('A', declared)}) == []


def test_suggest_sequence_finds_the_variant_that_matches():
    struc = read(TWO_CHAIN_PDB)
    candidates = {
        'default.fa': seq_dict(struc, 'A'),  # one chain only, so does not match
        '2-chain.fa': seq_dict(struc, 'A', 'B'),
    }
    assert sequence_check.check_sequences(struc, candidates['default.fa'])
    assert sequence_check.suggest_sequence(struc, candidates) == '2-chain.fa'


def test_suggest_sequence_returns_none_when_nothing_matches():
    struc = read(TWO_CHAIN_PDB)
    candidates = {'default.fa': {'A': ('A', 'ACDEFGHIK')}}
    assert sequence_check.suggest_sequence(struc, candidates) is None


def test_polymer_chains_ignores_waters_and_ligands():
    struc = read(ONE_CHAIN_PDB)
    # the file has a W chain of waters, which must not be treated as needing a sequence
    assert sorted(sequence_check.polymer_chains(struc)) == ['A']


def test_collect_candidates():
    sequences = {
        'default': 'default.fa',
        'variants': [{'sequence': '2-chain.fa', 'crystals': ['xtal-1', 'xtal-2']}],
    }
    default_seq = {'A': ('A', 'ACDEF')}
    variant_seq = {'A': ('A', 'ACDEF'), 'B': ('A', 'ACDEF')}
    candidates = sequence_check.collect_candidates(sequences, default_seq, {'xtal-1': variant_seq})
    assert candidates == {'default.fa': default_seq, '2-chain.fa': variant_seq}


def test_collect_candidates_without_sequences():
    assert sequence_check.collect_candidates(None, None, {}) == {}


def test_read_structure_reads_mmcif():
    struc = sequence_check.read_structure(mmcif_file='test-data/cif/disconnected_pair.cif')
    chains = sequence_check.polymer_chains(struc)
    assert sorted(chains) == ['A']
    assert sequence_check.check_sequences(struc, {name: ('A', observed(struc, name)) for name in chains}) == []


def test_structure_with_no_model_is_not_an_error():
    # a CIF that is not a set of coordinates at all, e.g. ligand restraints, must not blow up
    struc = sequence_check.read_structure(mmcif_file='test-data/inputs_2/ligand_bound_manual/Mpro-i0130.cif')
    assert sequence_check.polymer_chains(struc) == {}
