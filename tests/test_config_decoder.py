from xchemalign.decoder import decoder


def test_with_missing_file():
    # Arrange
    expected_error = "The config file 'no-such-file.yaml' does not exist"

    # Act
    error = decoder.validate_config_schema("no-such-file.yaml")

    # Assert
    assert error == expected_error


def test_test_data_config_1():
    # Arrange
    expected_error = None

    # Act
    error = decoder.validate_config_schema("test-data/config_1.yaml")

    # Assert
    assert error == expected_error


def test_test_data_config_2():
    # Arrange
    expected_error = None

    # Act
    error = decoder.validate_config_schema("test-data/config_2.yaml")

    # Assert
    assert error == expected_error


def test_test_data_config_3():
   # Arrange
    expected_error = None

    # Act
    error = decoder.validate_config_schema("test-data/config_3.yaml")

    # Assert
    assert error == expected_error
