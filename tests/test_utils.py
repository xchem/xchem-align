from xchemalign import utils


def test_parse_compound_smiles():
    inputs = {'string1': [1, 1], 'string1;string2': [2, 1, 1], 'string1;string21 string22;string3': [3, 1, 2, 1]}
    for s, r in inputs.items():
        result = utils.parse_compound_smiles(s)
        assert len(result) == r[0]
        for i, v in enumerate(result):
            assert len(v) == r[i + 1]
    print('OK')
