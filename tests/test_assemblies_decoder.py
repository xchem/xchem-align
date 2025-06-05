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


def test_test_data_assemblies_with_bad_assembly_name():
    # Arrange
    expected_error = "'dimer+1' does not match any of the regexes: '^[a-zA-Z0-9-_.]{1,80}$'"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-bad-assembly-name.yaml")

    # Assert
    assert error == expected_error


def test_test_data_assemblies_with_duplicate_crystalforms():
    # Arrange
    expected_error = "Found duplicate key 'xtalform1'"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-duplicate-crystalforms.yaml")

    # Assert
    assert error == expected_error


def test_test_data_assemblies_with_bad_assemblies_indentation():
    # Arrange
    expected_error = "Additional properties are not allowed ('1' was unexpected)"

    # Act
    error = decoder.validate_assemblies_schema("test-data/assemblies-with-bad-assemblies-indentation.yaml")

    # Assert
    assert error == expected_error
