DEBUG: bool = True

ALIGNABILITY_GRAPH_FILE_NAME: str = "alignability.gml"
CONNECTED_COMPONENTS_FILE_NAME: str = "connected_components.json"
TRANSFORMS_FILE_NAME: str = "transforms.json"
NEIGHBOURHOODS_FILE_NAME: str = "neighbourhoods.json"
DATA_JSON_PATH: str = "data.json"
SITES_TRANSFORMS_FILE_NAME: str = "site_transforms.json"
XTALFORMS_FILE_NAME: str = "xtalforms.json"
ASSIGNED_XTALFORMS_FILE_NAME: str = "assigned_xtalforms.json"
ASSEMBLIES_FILE_NAME: str = "assemblies.json"
CONFORMER_SITE_FILE: str = "conformer_sites.json"
CANONICAL_SITE_FILE: str = "canonical_sites.json"
XTALFORM_SITE_FILE: str = "ctalform_sites.json"

HIERARCHY_YAML = 'hierarchy.yaml'
BIOCHAIN_PRIORITIES_YAML = 'biochain_priorities.yaml'
ASSEMBLY_LANDMARKS_YAML = 'assembly_landmarks.yaml'
ASSEMBLY_TRANSFORMS_YAML = 'assembly_transforms.yaml'
CHAIN_TO_ASSEMBLY_YAML = 'chain_to_assembly.yaml'

OBSERVED_SITES_FILE_NAME: str = ""
SITES_FILE_NAME: str = "sites.json"

PANDDA_ANALYSES_DIR: str = "analyses"
PANDDA_EVENTS_INSPECT_TABLE_PATH: str = "pandda_inspect_events.csv"
PANDDA_FINAL_STRUCTURE_PDB_DIR: str = "modelled_structures"
PANDDA_FINAL_STRUCTURE_PDB_TEMPLATE: str = "{dtag}-pandda-model.pdb"

PANDDA_EVENT_MAP_TEMPLATE: str = "{dtag}-event_{event_id}_1-BDC_{bdc}_map.native.ccp4"


PANDDA_PROCESSED_DATASETS_DIR: str = "processed_datasets"

ALIGNED_STRUCTURES_DIR: str = "aligned_structures"
ALIGNED_XMAPS_DIR: str = "aligned_xmaps"
ALIGNED_FILES_DIR: str = "aligned_files"

MODEL_DIR_PDB: str = "refine.pdb"
MODEL_DIR_XMAP: str = "refine.ccp4"
MODEL_DIR_MTZ: str = "refine.mtz"

OUTPUT_JSON_PATH: str = "output.json"
ALIGNED_STRUCTURE_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}.pdb"
ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_artefacts.pdb"

ALIGNED_XMAP_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_sigmaa.ccp4"
ALIGNED_DIFF_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_diff.ccp4"
ALIGNED_EVENT_MAP_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_event.ccp4"

ALIGNED_XMAP_CRYSTALLOGRAPHIC_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_sigmaa_crystallographic.ccp4"
ALIGNED_DIFF_CRYSTALLOGRAPHIC_TEMPLATE: str = "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_diff_crystallographic.ccp4"
ALIGNED_EVENT_MAP_CRYSTALLOGRAPHIC_TEMPLATE: str = (
    "{dtag}_{chain}_{residue}_{altloc}_{version}_{site}_event_crystallographic.ccp4"
)


FS_MODEL_YAML_FILE_NAME = "fs_model.yaml"
ASSEMBLIES_YAML_FILE_NAME = "assemblies.yaml"
XTALFORMS_YAML_FILE_NAME = "xtalforms.yaml"
ASSIGNED_XTALFORMS_YAML_FILE_NAME = "assigned_xtalforms.yaml"
NEIGHBOURHOODS_YAML_FILE_NAME = "neighbourhoods.yaml"
CONNECTED_COMPONENTS_YAML_NAME = "connected_components.yaml"
TRANSFORMS_YAML_FILE_NAME = "neighbourhood_transforms.yaml"
CONFORMER_SITE_YAML_FILE = "conformer_sites.yaml"
CONFORMER_SITES_TRANSFORMS_YAML_FILE_NAME = "conformer_site_transforms.yaml"
CANONICAL_SITE_YAML_FILE = "canonical_sites.yaml"
CANONICAL_SITES_TRANSFORMS_YAML_FILE_NAME = "canonical_site_transforms.yaml"
XTALFORM_SITE_YAML_FILE = "xtalform_sites.yaml"
REFERENCE_STRUCTURE_TRANSFORMS_YAML = "reference_structure_transforms.yaml"

ALIGNED_REFERENCE_STRUCTURE_TEMPLATE = "{dtag}_{site}_ref.pdb"
ALIGNED_REFERENCE_STRUCTURE_ARTEFACTS_TEMPLATE = "{dtag}_{site}_artefacts_ref.pdb"
ALIGNED_REFERENCE_XMAP_TEMPLATE = "{dtag}_{site}_ref.ccp4"

RESIDUE_NAMES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "GLY",
]
