"""A module to validate and decode XChemAlign definitions.

This decoder module offers validation of 'assemblies' and 'config' files
using YAML BaseLoader and a JSON schema to enforce the structure of the file.
Following these checks there are additional (custom) 'content' tests to validate parts
of the file that are not possible using a YAML loader or JSON schema.

Assembly validation is done by calling 'validate_assemblies_schema()' with an assemblies
file. This returns an error string if there's an error or None if the file is OK.

Config validation is done by calling 'validate_configs_schema()' with a config
file. This returns an error string if there's an error or None if the file is OK.

Things that are checked (Assemblies):

1.  Structure - handled by YAML BaseLoader
    a.  Detects duplicate keys.
        e.g. 'xtalform1' being used twice in the crystalform block

2.  Structure - handled by the jsonschema validation
    a.  Does it contain all the expected keys in the right place?
        like 'assemblies', 'biomol', and 'chains'
    b.  Does it contain any unexpected keys?
    c.  Is the assembly key valid? (up to 80 characters using [a-zA-Z0-9-_.])
        e.g. 'dimer'
    d.  Is the crystalform key valid?
        e.g. 'xtalform3'

3.  Content - handled by the custom function _validate_assemblies_content()
    a.  Does every crystalform assembly refer to an assembly in the file?

Things that are checked (Config):

1.  Structure - handled by YAML BaseLoader
    a.  Detects duplicate keys.

2.  Structure - handled by the jsonschema validation
    a.  Does it contain all the expected keys in the right place?
    b.  Does it contain any unexpected keys?

3.  Content - handled by the custom function _validate_config_content()
    a.  Nothing additional is checked atm

Ideally user should add a reference to the schema in their YAML files.
If they do this the editor (if sufficiently empowered) should provide live
feedback on basic YAML formatting deviations. The following at the start of the
file, as an example (shortened for clarity), in VisualStudioCode will
allow the editor to display live errors with Assemblies file for example: -

# yaml-language-server: $schema=https://raw.githubusercontent.com/xchem/xchem-align/[...]/assemblies-schema.yaml
"""

import os
from typing import Any

import jsonschema
import yaml
from yaml.constructor import ConstructorError

# The (built-in) schemas...
# from the same directory as us.
_ASSEMBLIES_SCHEMA_FILE: str = os.path.join(os.path.dirname(__file__), "assemblies-schema.yaml")
_CONFIG_SCHEMA_FILE: str = os.path.join(os.path.dirname(__file__), "config-schema.yaml")

# Load the Assemblies schema YAML file now.
# This must work as the file is installed along with this module.
assert os.path.isfile(_ASSEMBLIES_SCHEMA_FILE)
with open(_ASSEMBLIES_SCHEMA_FILE, "r", encoding="utf8") as schema_file:
    _ASSEMBLIES_SCHEMA: dict[str, Any] = yaml.load(schema_file, Loader=yaml.FullLoader)
assert _ASSEMBLIES_SCHEMA

# Load the Config schema YAML file now.
# This must work as the file is installed along with this module.
assert os.path.isfile(_CONFIG_SCHEMA_FILE)
with open(_CONFIG_SCHEMA_FILE, "r", encoding="utf8") as schema_file:
    _CONFIG_SCHEMA: dict[str, Any] = yaml.load(schema_file, Loader=yaml.FullLoader)
assert _CONFIG_SCHEMA


# A YAML constructor and custom BaseLoader class
# that detects duplicate YAML keys
def _NO_DUPLICATES_CONSTRUCTOR(loader, node, deep=False):
    """Check for duplicate keys."""
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        value = loader.construct_object(value_node, deep=deep)
        if key in mapping:
            raise ConstructorError(f"Found duplicate key '{key}'")
        mapping[key] = value

    return loader.construct_mapping(node, deep)


class DupCheckLoader(yaml.BaseLoader):
    """Local class to prevent pollution of global yaml.Loader."""

    pass


DupCheckLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _NO_DUPLICATES_CONSTRUCTOR)


def validate_assemblies_schema(assembly_filename: str) -> str | None:
    """Checks the Assemblies definition against the built-in schema.
    If there's an error the error text is returned, otherwise None.
    """
    assert assembly_filename
    assemblies, error = _validate_file_schema(input_filename=assembly_filename, schema=_ASSEMBLIES_SCHEMA)
    return error or _validate_assemblies_content(assemblies)


def validate_config_schema(config_filename: str) -> str | None:
    """Checks the Config definition against the built-in schema.
    If there's an error the error text is returned, otherwise None.
    """
    assert config_filename
    config, error = _validate_file_schema(input_filename=config_filename, schema=_CONFIG_SCHEMA)
    return error or _validate_config_content(config)


def _validate_file_schema(*, input_filename: str, schema: dict[str, Any]) -> tuple[dict[str, Any] | None, str | None]:
    """Does schema validation for a file, given a schema."""

    if not os.path.isfile(input_filename):
        return None, "The file does not exist"

    try:
        with open(input_filename, "r", encoding="utf8") as yaml_file:
            input_yaml: dict[str, Any] = yaml.load(yaml_file, DupCheckLoader)
    except yaml.constructor.ConstructorError as cex:
        return None, str(cex)
    except yaml.scanner.ScannerError:
        return None, "Unable to understand the file - content is not valid YAML"

    if not input_yaml:
        return None, "The file is empty"

    try:
        jsonschema.validate(input_yaml, schema=schema)
    except jsonschema.ValidationError as vex:
        return None, str(vex.message)
    except TypeError as tex:
        return None, str(tex)

    # OK if we get here
    # Return the loaded YAML map.
    return input_yaml, None


def _validate_assemblies_content(assemblies_content: dict[str, Any]) -> str | None:
    """Assuming the file has already passed schema validation this function
    checks additional content, like cross-references of assemblies."""
    assert isinstance(assemblies_content, dict)
    # Check assembly cross-references
    # What assemblies are declared?
    assemblies: list[str] = []
    assemblies.extend(iter(assemblies_content["assemblies"]))
    # Does each crystalform refer to one of the known assemblies?
    crystalforms = assemblies_content["crystalforms"]
    for crystalform in crystalforms:
        crystalfrom_assemblies = crystalforms[f"{crystalform}"]["assemblies"]
        for assembly in crystalfrom_assemblies:
            assembly_name: str = crystalfrom_assemblies[f"{assembly}"]["assembly"]
            if assembly_name not in assemblies:
                return f"The assembly '{assembly_name}' in crystalform '{crystalform}->{assembly}' is not an assembly in the file"

    # OK if we get here
    return None


def _validate_config_content(config_content: dict[str, Any]) -> str | None:
    """Assuming the file has already passed schema validation this function
    checks additional content."""
    assert isinstance(config_content, dict)

    # OK if we get here
    return None
