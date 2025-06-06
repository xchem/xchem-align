import glob

from xchemalign.decoder import decoder

_SKIP: list[str] = []


def test_all_test_data_assemblies_files():
    # Search for all the assemblies files (in test-data) and test them
    num_files_tested = 0
    config_files = glob.glob("test-data/**/assemblies.yaml", recursive=True)
    for config_file in config_files:
        if config_file not in _SKIP:
            error = decoder.validate_assemblies_schema(config_file)
            num_files_tested += 1
            assert error is None, config_file
    assert num_files_tested > 0
