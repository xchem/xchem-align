from pathlib import Path

from ligand_neighbourhood_alignment.data import Transforms


def save_transforms(transforms: Transforms, path: Path):
    with open(path, "w") as f:
        f.write(transforms.json())
