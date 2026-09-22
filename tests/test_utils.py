from pathlib import Path

from rdkit import Chem

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


def test_sequence_file_paths_defaults():
    # sequences.dir and sequences.default both omitted, so the documented defaults apply
    input_yaml = {
        Constants.CONFIG_DIR: 'some/dir',
        Constants.CONFIG_SEQUENCES: {
            Constants.CONFIG_VARIANTS: [
                {Constants.CONFIG_SEQUENCE: '2-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-1']}
            ]
        },
    }
    assert [str(p) for p in utils.sequence_file_paths(input_yaml)] == [
        'processing/analysis/sequences/default.fa',
        'processing/analysis/sequences/2-chain.fa',
    ]


def test_sequence_file_paths_with_an_empty_sequences_section():
    # an empty section means the same as no section, matching read_sequences
    input_yaml = {Constants.CONFIG_DIR: 'some/dir', Constants.CONFIG_SEQUENCES: {}}
    assert utils.sequence_file_paths(input_yaml) == []
    assert utils.read_sequences(Path('/nonexistent'), input_yaml, fatal=False) == (None, {})


def test_sequence_file_paths_with_variants():
    input_yaml = {
        Constants.CONFIG_DIR: 'some/dir',
        Constants.CONFIG_SEQUENCES: {
            Constants.CONFIG_DIR: 'processing/sequences',
            Constants.CONFIG_DEFAULT: 'default.fa',
            Constants.CONFIG_VARIANTS: [
                {Constants.CONFIG_SEQUENCE: '2-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-1']},
                {Constants.CONFIG_SEQUENCE: '3-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-2']},
            ],
        },
    }
    assert [str(p) for p in utils.sequence_file_paths(input_yaml)] == [
        'processing/sequences/default.fa',
        'processing/sequences/2-chain.fa',
        'processing/sequences/3-chain.fa',
    ]


def test_sequence_file_paths_deduplicates():
    # the same file may legitimately be named by more than one variant
    input_yaml = {
        Constants.CONFIG_SEQUENCES: {
            Constants.CONFIG_VARIANTS: [
                {Constants.CONFIG_SEQUENCE: '2-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-1']},
                {Constants.CONFIG_SEQUENCE: '2-chain.fa', Constants.CONFIG_CRYSTALS: ['xtal-2']},
            ]
        }
    }
    paths = [str(p) for p in utils.sequence_file_paths(input_yaml)]
    assert paths == ['processing/analysis/sequences/default.fa', 'processing/analysis/sequences/2-chain.fa']


def test_sequence_file_paths_without_a_sequences_section():
    assert utils.sequence_file_paths({Constants.CONFIG_DIR: 'some/dir'}) == []


def test_parse_compound_smiles():
    inputs = {'string1': [1, 1], 'string1;string2': [2, 1, 1], 'string1;string21 string22;string3': [3, 1, 2, 1]}
    for s, r in inputs.items():
        result = utils.parse_compound_smiles(s)
        assert len(result) == r[0]
        for i, v in enumerate(result):
            assert len(v) == r[i + 1]
    print('OK')


def write_ligand_cif(tmp_path, charge_column, atoms, bonds):
    """
    Write a single-ligand CIF without hydrogens (RDKit supplies them implicitly).

    :param charge_column: 'charge' or 'partial_charge'
    :param atoms: list of (atom_id, element, charge) tuples
    :param bonds: list of (atom_id_1, atom_id_2, type) tuples
    """
    lines = [
        'data_comp_LIG',
        'loop_',
        '_chem_comp_atom.comp_id',
        '_chem_comp_atom.atom_id',
        '_chem_comp_atom.type_symbol',
        '_chem_comp_atom.' + charge_column,
        '_chem_comp_atom.x',
        '_chem_comp_atom.y',
        '_chem_comp_atom.z',
    ]
    for i, (atom_id, element, charge) in enumerate(atoms):
        lines.append('LIG {} {} {} {:.3f} {:.3f} 0.000'.format(atom_id, element, charge, 1.5 * i, 0.3 * (i % 2)))
    lines += ['loop_', '_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1', '_chem_comp_bond.atom_id_2']
    lines.append('_chem_comp_bond.type')
    for a1, a2, bond_type in bonds:
        lines.append('LIG {} {} {}'.format(a1, a2, bond_type))
    p = tmp_path / 'LIG.cif'
    p.write_text('\n'.join(lines) + '\n')
    return p


def test_gen_mols_from_cif_whole_charges_are_kept(tmp_path):
    p = write_ligand_cif(tmp_path, 'charge', [('C1', 'C', '0'), ('N1', 'N', '1')], [('C1', 'N1', 'single')])
    mol = utils.gen_mols_from_cif(p)[0]
    assert Chem.MolToSmiles(mol) == 'C[NH3+]'


def test_gen_mols_from_cif_charge_spread_over_deloc_group_is_kept(tmp_path):
    # old CCP4 monomer-library style carboxylate: -0.5 on each O, which round() turned into 0
    p = write_ligand_cif(
        tmp_path,
        'partial_charge',
        [('C2', 'C', '0.000'), ('C1', 'C', '0.000'), ('O1', 'O', '-0.500'), ('O2', 'O', '-0.500')],
        [('C2', 'C1', 'single'), ('C1', 'O1', 'deloc'), ('C1', 'O2', 'deloc')],
    )
    mol = utils.gen_mols_from_cif(p)[0]
    assert Chem.GetFormalCharge(mol) == -1
    charges = {a.GetProp('atom_id'): a.GetFormalCharge() for a in mol.GetAtoms()}
    assert charges == {'C2': 0, 'C1': 0, 'O1': -1, 'O2': 0}
    # deloc is still read as single (issue #103), so the other O is not yet a double bond
    assert Chem.MolToSmiles(mol) == 'CC([O-])O'


def test_gen_mols_from_cif_partial_charges_are_not_formal_charges(tmp_path):
    # rounding used to turn these into +1 on the carbonyl C and -1 on the O, which failed sanitization
    p = write_ligand_cif(
        tmp_path,
        'partial_charge',
        [('C1', 'C', '-0.100'), ('C2', 'C', '0.600'), ('O1', 'O', '-0.600'), ('C3', 'C', '-0.100')],
        [('C1', 'C2', 'single'), ('C2', 'O1', 'double'), ('C2', 'C3', 'single')],
    )
    mol = utils.gen_mols_from_cif(p)[0]
    assert Chem.MolToSmiles(mol) == 'CC(C)=O'
