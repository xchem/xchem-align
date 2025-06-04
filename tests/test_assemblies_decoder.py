from xchemalign.decoder import decoder


def test_test_data_assemblies():
    # Arrange
    expected_error = None

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies.yaml")

    # Assert
    assert error == expected_error


def test_test_data_lna_assemblies():
    # Arrange
    expected_error = None

    # Act
    error = decoder.validate_assemblies_schema("test-data/lna_assemblies.yaml")

    # Assert
    assert error == expected_error


def test_test_data_assemblies_with_missing_crystalform_assembly():
    # Arrange
    expected_error = "The assembly 'missing-dimer' in crystalform 'xtalform2->1' is not an assembly in the file"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-missing-crystalform-assembly.yaml")

    # Assert
    assert error == expected_error
