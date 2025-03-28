from pathlib import Path

from ligand_neighbourhood_alignment.data import LigandNeighbourhoods


def save_neighbourhoods(ligand_neighbourhoods: LigandNeighbourhoods, path: Path):
    with open(path, "w") as f:
        f.write(ligand_neighbourhoods.json())
