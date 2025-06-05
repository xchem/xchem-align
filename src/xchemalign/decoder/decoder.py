"""A module to validate and decode XChemAlign definitions.

This decoder module offers validation of 'assemblies' files using YAML BaseLoader
and a JSON schema to enforce the structure of the file. Following these checks
there are additional (custom) 'content' tests to validate parts of the file
that are not possible using a YAML loader or JSON schema.

Assembly validation is done by calling 'validate_assemblies_schema()' with an assemblies
file. This returns an error string if there's an error or None if the file is OK.

Things that are checked:

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

Ideally user should add a reference to the schema in their YAML files.
If they do this the editor (if sufficiently empowered) should provide live
feedback on basic YAML formatting deviations. The following at the start of the
file, as an example (shortened for clarity), in VisualStudioCode will
allow the editor to display live errors: -

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

# Load the Workflow schema YAML file now.
# This must work as the file is installed along with this module.
assert os.path.isfile(_ASSEMBLIES_SCHEMA_FILE)
with open(_ASSEMBLIES_SCHEMA_FILE, "r", encoding="utf8") as schema_file:
    _ASSEMBLIES_SCHEMA: dict[str, Any] = yaml.load(schema_file, Loader=yaml.FullLoader)
assert _ASSEMBLIES_SCHEMA


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
    if not os.path.isfile(assembly_filename):
        return f"The assembly file '{assembly_filename}' does not exist"

    try:
        with open(assembly_filename, "r", encoding="utf8") as assembly_file:
            assembly: dict[str, Any] = yaml.load(assembly_file, DupCheckLoader)
    except yaml.constructor.ConstructorError as cex:
        return str(cex)

    try:
        jsonschema.validate(assembly, schema=_ASSEMBLIES_SCHEMA)
    except jsonschema.ValidationError as vex:
        return str(vex.message)
    except TypeError as tex:
        return str(tex)

    # OK so far, now check additional content,
    # i.e. does each crystalform assembly refer to an assembly?
    return _validate_assemblies_content(assembly)


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
