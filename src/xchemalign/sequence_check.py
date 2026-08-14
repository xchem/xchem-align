"""
Checks that the amino acid sequences a user declares in config.yaml actually match the residues
present in the refined models.

A wrong sequence is not rejected until the structures reach the wwPDB, so it is worth catching
early. See https://github.com/m2ms/fragalysis-frontend/issues/2342

The residue comparison has to be an alignment rather than a string comparison because:
  - model residue numbering does not correspond to position in the declared sequence (models
    routinely start at residue 7, 486, 582 ...);
  - residues that were not modelled leave gaps, which is normal and must not be reported;
  - a model can contain an inserted residue, which shifts everything after it.

gemmi already bundles an aligner, so no extra dependency is needed for this.
"""

import gemmi

from xchemalign.utils import Constants

# 's' has gentle, symmetric penalties. Not 'p' ('partial model'): its mismatch and gap-open costs of
# -10000 rely on being given a per-residue gap-open array derived from the numbering, and without one
# a single mismatched residue collapses the whole alignment.
SCORING = gemmi.AlignmentScoring('s')

# what expand_one_letter_sequence/find_tabulated_residue use for anything it cannot identify
UNKNOWN_RESIDUE = 'X'


def polymer_chains(struc):
    """
    The peptide chains of a structure, keyed by chain name. Matches how adjust_chains_and_entities in
    pdb_deposition decides what is a protein chain, so the two agree on which chains need a sequence.

    :param struc: the gemmi Structure
    :return: dict of chain name -> gemmi ResidueSpan
    """
    chains: dict = {}
    if len(struc) == 0:
        return chains
    for chain in struc[0]:  # always only a single model
        polymer = chain.get_polymer()
        if len(polymer) > 0 and polymer.check_polymer_type() == gemmi.PolymerType.PeptideL:
            chains[chain.name] = polymer
    return chains


def observed_sequence(polymer):
    """
    The one letter sequence of the residues actually present, with modified residues represented by
    the residue they derive from, so that e.g. selenomethionine counts as the M the user declared.

    :param polymer: a gemmi ResidueSpan
    :return: one letter sequence as a string
    """
    codes = []
    for residue in polymer:
        info = gemmi.find_tabulated_residue(residue.name)
        # one_letter_code is lower case for non-standard residues that have a standard parent
        code = info.one_letter_code.upper() if info else ''
        codes.append(code if code.isalpha() else UNKNOWN_RESIDUE)
    return ''.join(codes)


def check_sequences(struc, seq_dict):
    """
    Check a structure against the sequences declared for it.

    :param struc: the gemmi Structure
    :param seq_dict: dict keyed by chain name, values a tuple of (entity name, one letter sequence),
        as returned by utils.read_fasta
    :return: list of descriptions of what does not match. Empty means the sequences are good.
    """
    issues = []
    chains = polymer_chains(struc)

    for name in sorted(set(seq_dict) - set(chains)):
        issues.append('chain {} has a declared sequence but is not in the model'.format(name))
    for name in sorted(set(chains) - set(seq_dict)):
        issues.append('chain {} is in the model but has no declared sequence'.format(name))

    for name in sorted(set(seq_dict) & set(chains)):
        declared = seq_dict[name][1]
        if not declared:
            issues.append('chain {} has an empty declared sequence'.format(name))
            continue
        issues.extend(compare_chain(name, declared, chains[name]))

    return issues


def compare_chain(name, declared, polymer):
    """
    Align a declared sequence to the residues of one chain and describe anything that does not fit.

    Residues that are in the declared sequence but not modelled are expected and are not reported;
    the check is that every residue that *is* modelled is accounted for by the declared sequence.

    :param name: the chain name, used in the descriptions
    :param declared: the declared one letter sequence
    :param polymer: a gemmi ResidueSpan
    :return: list of descriptions
    """
    issues = []
    observed = observed_sequence(polymer)
    alignment = gemmi.align_string_sequences(list(declared), list(observed), [], SCORING)

    declared_aligned = alignment.add_gaps(declared, 1)
    observed_aligned = alignment.add_gaps(observed, 2)

    index = 0  # index into polymer, advanced for each aligned position that is a real residue
    for declared_code, observed_code in zip(declared_aligned, observed_aligned):
        if observed_code == '-':
            # a residue that was not modelled, which is normal
            continue
        residue = polymer[index]
        index += 1
        if declared_code == '-':
            issues.append(
                'chain {} has {} {} in the model, which the declared sequence does not account for'.format(
                    name, residue.name, residue.seqid.num
                )
            )
        elif declared_code != observed_code:
            issues.append(
                'chain {} declares {} where the model has {} {}'.format(
                    name, declared_code, residue.name, residue.seqid.num
                )
            )

    return issues


def collect_candidates(sequences_config, default_seq, variants):
    """
    The sequences declared for an input, keyed by the FASTA file they came from, so that a crystal that
    fails its own sequence can be told which of the others it does match.

    :param sequences_config: the sequences section of the input in config.yaml
    :param default_seq: the default chain->(entity, seq) dict, as returned by utils.read_sequences
    :param variants: the crystal name -> variant dict, as returned by utils.read_sequences
    :return: dict of FASTA file name -> seq_dict
    """
    candidates: dict = {}
    if not sequences_config or default_seq is None:
        return candidates

    candidates[sequences_config.get(Constants.CONFIG_DEFAULT, 'default.fa')] = default_seq
    for variant in sequences_config.get(Constants.CONFIG_VARIANTS) or []:
        name = variant.get(Constants.CONFIG_SEQUENCE)
        for xtal in variant.get(Constants.CONFIG_CRYSTALS) or []:
            # every crystal of a variant maps to the same dict, so the first one that was read is enough
            if xtal in variants:
                candidates[name] = variants[xtal]
                break
    return candidates


def suggest_sequence(struc, candidates):
    """
    Find a declared sequence that would have matched this structure. Used to turn 'this crystal is
    wrong' into 'this crystal matches 2-chain.fa', which is the fix the variants section exists for.

    :param struc: the gemmi Structure
    :param candidates: dict of name (the FASTA file name) -> seq_dict
    :return: the name of one that matches, or None
    """
    for name in sorted(candidates):
        seq_dict = candidates[name]
        if seq_dict and not check_sequences(struc, seq_dict):
            return name
    return None


def read_structure(pdb_file=None, mmcif_file=None):
    """
    Read a model the same way the deposition pipeline does, so the checks see the same residues.

    :param pdb_file: path of a PDB file, used if no mmcif_file is given
    :param mmcif_file: path of an mmCIF file, preferred if present
    :return: the gemmi Structure
    """
    if mmcif_file:
        struc = gemmi.read_structure(str(mmcif_file))
    else:
        struc = gemmi.read_pdb(str(pdb_file), ignore_ter=True)
    struc.setup_entities()
    return struc
