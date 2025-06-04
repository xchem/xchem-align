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
    with open(assembly_filename, "r", encoding="utf8") as assembly_file:
        assembly: dict[str, Any] = yaml.load(assembly_file, Loader=yaml.BaseLoader)

    try:
        jsonschema.validate(assembly, schema=_ASSEMBLIES_SCHEMA)
    except jsonschema.ValidationError as vex:
        return str(vex.message)
    except TypeError as tex:
        return str(tex)

    # OK if we get here
    return None
