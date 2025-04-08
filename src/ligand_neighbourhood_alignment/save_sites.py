from pathlib import Path

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.data import CanonicalSites


def save_sites(sites: CanonicalSites, path: Path):
    with open(path / constants.SITES_FILE_NAME, "w") as f:
        f.write(sites.json())
