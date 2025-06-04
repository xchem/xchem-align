from xchemalign.decoder import decoder


def test_with_missigg_file():
    # Arrange
    expected_error = "The assembly file 'no-such-file.yaml' does not exist"

    # Act
    error = decoder.validate_assemblies_schema("no-such-file.yaml")

    # Assert
    assert error == expected_error


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


def test_test_data_assemblies_with_missing_assemblies_section():
    # Arrange
    expected_error = "'assemblies' is a required property"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-missing-assemblies-section.yaml")

    # Assert
    assert error == expected_error


def test_test_data_assemblies_with_mis_spelled_assemblies_section():
    # Arrange
    expected_error = "Additional properties are not allowed ('azzemblies' was unexpected)"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-mis-spelled-assemblies-section.yaml")

    # Assert
    assert error == expected_error
