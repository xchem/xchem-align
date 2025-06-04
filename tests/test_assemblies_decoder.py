from xchemalign.decoder import decoder

def test_test_data_assemblies():
    # Arrange

    # Act
    error = decoder.validate_assemblies_schema('test-data/assemblies.yaml')

    # Assert
    assert error is None


def test_test_data_lna_assemblies():
    # Arrange

    # Act
    error = decoder.validate_assemblies_schema('test-data/lna_assemblies.yaml')

    # Assert
    assert error is None
