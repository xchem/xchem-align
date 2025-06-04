"""A module to validate and decode XChemAlign definitions.
"""

import os
from typing import Any

import jsonschema
import yaml

# The (built-in) schemas...
# from the same directory as us.
_ASSEMBLIES_SCHEMA_FILE: str = os.path.join(os.path.dirname(__file__), "assemblies-schema.yaml")

# Load the Workflow schema YAML file now.
# This must work as the file is installed along with this module.
assert os.path.isfile(_ASSEMBLIES_SCHEMA_FILE)
with open(_ASSEMBLIES_SCHEMA_FILE, "r", encoding="utf8") as schema_file:
    _ASSEMBLIES_SCHEMA: dict[str, Any] = yaml.load(schema_file, Loader=yaml.FullLoader)
assert _ASSEMBLIES_SCHEMA


def validate_assemblies_schema(assembly_filename: str) -> str | None:
    """Checks the Assemblies definition against the built-in schema.
    If there's an error the error text is returned, otherwise None.
    """
    assert assembly_filename
    if not os.path.isfile(assembly_filename):
        return f"The assembly file does not exist ({assembly_filename})"

    with open(assembly_filename, "r", encoding="utf8") as assembly_file:
        assembly: dict[str, Any] = yaml.load(assembly_file, Loader=yaml.BaseLoader)

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
    assemblies.extend(iter(assemblies_content['assemblies']))
    # Does each crystalform refer to one of the known assemblies?
    crystalforms = assemblies_content['crystalforms']
    for crystalform in crystalforms:
        crystalfrom_assemblies = crystalforms[f'{crystalform}']['assemblies']
        for assembly in crystalfrom_assemblies:
            assembly_name: str = crystalfrom_assemblies[f'{assembly}']['assembly']
            if assembly_name not in assemblies:
                return f"The assembly '{assembly_name}' in crystalform '{crystalform}->{assembly}' is not an assembly in this file"

    # OK if we get here
    return None
