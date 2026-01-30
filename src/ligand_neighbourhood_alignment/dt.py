import json
import os
import re
from pathlib import Path
from typing import Generator
import sys

import pandas as pd
import yaml
from loguru import logger
import gemmi
import networkx as nx
from loguru import logger
from pydantic import BaseModel, field_validator

from ligand_neighbourhood_alignment import constants



logger.remove()  # for someone not familiar with the lib, whats going on here?
logger.add(sys.stdout, level="INFO")

AlignmentHeirarchy = dict[str, tuple[str, str]]
StructureLandmarks = dict[tuple[str, tuple[str, str], str], tuple[float, float, float]]
Structure = gemmi.Structure
LigandNeighbourhoodID = tuple[str,str,str,str,str]

class DatasetID(BaseModel):
    dtag: str

    def __hash__(self):
        return hash(self.dtag)

    def __eq__(self, other):
        try:
            if other.dtag == self.dtag:
                return True
        except Exception:
            if self.dtag == other:
                return True
            else:
                return False


class LigandID(BaseModel):
    dtag: str
    chain: str
    residue: int

    def __eq__(self, other) -> bool:
        try:
            if self.dtag == other.dtag:
                if self.residue == other.residue:
                    if self.chain == other.chain:
                        return True
            return False
        except Exception:
            return False

    def __hash__(self):
        return hash((self.dtag, self.chain, self.residue))

    def to_string(
        self,
    ):
        return f"{self.dtag}~{self.chain}~{self.residue}"

    @classmethod
    def from_string(cls, string):
        dtag, chain, residue = string.split("~")

        return LigandID(dtag=dtag, chain=chain, residue=int(residue))


class SymOp(BaseModel):
    operation: str
    image: tuple[int, int, int]


class AtomID(BaseModel):
    chain: str
    residue: int
    atom: str

    def __eq__(self, other) -> bool:
        if self.atom == other.atom:
            if self.residue == other.residue:
                if self.chain == other.chain:
                    return True

        return False

    def __hash__(self):
        return hash((self.chain, self.residue, self.atom))


class Transform(BaseModel):
    vec: list[float]
    mat: list[list[float]]

    def __hash__(self):
        return hash(str([x for x in self.vec] + [x for y in self.mat for x in y]))


class Transforms(BaseModel):
    ligand_ids: list[tuple[LigandID, LigandID]]
    transforms: list[Transform]

    def get_transform(self, transform_id: tuple[LigandID, LigandID]):
        for ligand_id, _transform in zip(self.ligand_ids, self.transforms):
            if transform_id == ligand_id:
                transform = gemmi.Transform()
                transform.vec.fromlist(_transform.vec)
                transform.mat.fromlist(_transform.mat)

                return transform

        raise Exception(f"Transform {transform_id} not in transforms {self.ligand_ids}!")


class Atom(BaseModel):
    element: str
    atom_id: AtomID
    x: float
    y: float
    z: float
    image: Transform


class AlignableSite(BaseModel):
    id: int
    name: str
    ligand_ids: list[LigandID]
    reference: LigandID


class AssemblyGenerator(BaseModel):
    id: int
    reference_chain: str
    chain: str
    triplet: str


class Assembly(BaseModel):
    id: int
    reference_assembly: int
    reference: DatasetID
    assembly: dict[int, AssemblyGenerator]


class Assemblies(BaseModel):
    assemblies: dict[int, Assembly]

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


class XtalFormAssembly(BaseModel):
    id: int
    reference_assembly: int
    generators: dict[int, AssemblyGenerator]


class XtalForm(BaseModel):
    id: int
    reference: DatasetID
    assemblies: dict[int, XtalFormAssembly]


class XtalForms(BaseModel):
    xtalforms: dict[int, XtalForm]

    def iter(self):
        for xtalform_id, xtalform in self.xtalforms.items():
            yield xtalform_id, xtalform

    def get_xtalform(self, item):
        for xtalform_id, xtalform in self.iter():
            if xtalform_id == item:
                return xtalform

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


class DatasetXtalforms(BaseModel):
    dataset_xtalforms: dict[str, int]


class XtalFormSite(BaseModel):
    id: int
    site_id: int
    xtalform_id: int
    crystallographic_chain: str
    members: list[LigandID]


class XtalFormSites(BaseModel):
    xtalform_sites: dict[int, XtalFormSite]

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


class SiteObservation(BaseModel):
    id: int
    ligand_id: LigandID
    xtal_form_site_id: int
    fragalysis_label: str
    description: str
    code: str
    dataset: str
    # compound: int


class LigandBindingEvent(BaseModel):
    id: int
    dtag: str
    chain: str
    residue: int
    xmap: str


class LigandBindingEvents(BaseModel):
    ligand_ids: list[LigandID]
    ligand_binding_events: list[LigandBindingEvent]

    def __getitem__(self, lid: LigandID) -> LigandBindingEvent:
        for _lid, lbe in zip(self.ligand_ids, self.ligand_binding_events):
            if lid == _lid:
                return lbe

        raise Exception(f"{lid} : {self.ligand_ids}")


class Dataset(BaseModel):
    dtag: str
    pdb: str
    xmap: str
    mtz: str
    ligand_binding_events: LigandBindingEvents


class Datasource(BaseModel):
    path: str
    datasource_type: str


class PanDDA(BaseModel):
    path: str
    event_table_path: str


class SystemData(BaseModel):
    datasources: list[Datasource]
    panddas: list[PanDDA]

    dataset_ids: list[DatasetID]
    datasets: list[Dataset]

    def get_dataset(self, did: str | DatasetID) -> Dataset:
        for _did, dataset in zip(self.dataset_ids, self.datasets):
            if _did == did:
                return dataset

        raise Exception(f"{did} : {self.dataset_ids}")

    def iter(self):
        for dataset_id, dataset in zip(self.dataset_ids, self.datasets):
            yield dataset_id, dataset


class LigandNeighbourhood(BaseModel):
    atom_ids: list[AtomID]
    atoms: list[Atom]
    artefact_atom_ids: list[AtomID]
    artefact_atoms: list[Atom]


class LigandNeighbourhoods(BaseModel):
    ligand_ids: list[LigandID]
    ligand_neighbourhoods: list[LigandNeighbourhood]

    def get_neighbourhood(self, ligand_id: LigandID):
        for _ligand_id, _neighbourhood in zip(self.ligand_ids, self.ligand_neighbourhoods):
            if _ligand_id == ligand_id:
                return _neighbourhood


# class Block(BaseModel):
#     xi: int
#     yi: int
#     zi: int
#     xmi: int
#     ymi: int
#     zmi: int
#     dx: int
#     dy: int
#     dz: int
#     transform: Transform

class Block:
    def __init__(
        self,
        xi: int,
        yi: int,
        zi: int,
        xmi: int,
        ymi: int,
        zmi: int,
        dx: int,
        dy: int,
        dz: int,
        transform: Transform
    ):
        self.xi: int = xi
        self.yi: int = yi
        self.zi: int = zi
        self.xmi: int = xmi
        self.ymi: int = ymi
        self.zmi: int = zmi
        self.dx: int = dx
        self.dy: int = dy
        self.dz: int = dz 
        self.transform: Transform = transform

class ResidueID(BaseModel):
    chain: str
    residue: int

    def __hash__(self):
        return hash((self.chain, self.residue))


class ConformerSite(BaseModel):
    id: int
    name: str
    residues: list[ResidueID]
    members: list[LigandID]
    reference_ligand_id: LigandID


class ConformerSites(BaseModel):
    conformer_sites: dict[int, ConformerSite]

    def iter(self) -> Generator[tuple[int, ConformerSite], None, None]:
        for cs_id, cs in self.conformer_sites.items():
            yield cs_id, cs

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


class CanonicalSite(BaseModel):
    id: int
    subsite_ids: list[int]
    subsites: list[ConformerSite]
    members: list[LigandID]
    residues: list[ResidueID]
    reference_ligand_id: LigandID
    reference_subsite_id: int
    reference_subsite: ConformerSite

    def iter(self) -> Generator[tuple[int, ConformerSite], None, None]:
        for subsite_id, subsite in zip(self.subsite_ids, self.subsites):
            yield subsite_id, subsite

    def get_subsite(self, subsite_id: int):
        for _subsite_id, subsite in self.iter():
            if subsite_id == _subsite_id:
                return subsite

        raise Exception(f"Site {subsite_id} not in sites: {self.subsite_ids}")


class CanonicalSites(BaseModel):
    site_ids: list[int]
    sites: list[CanonicalSite]
    reference_site: CanonicalSite
    reference_site_id: int

    def iter(self) -> Generator[tuple[int, CanonicalSite], None, None]:
        for site_id, site in zip(self.site_ids, self.sites):
            yield site_id, site

    def get_site(self, site_id: int):
        for _site_id, site in self.iter():
            if site_id == _site_id:
                return site

        raise Exception(f"Site {site_id} not in sites: {self.site_ids}")

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


class SystemSites(BaseModel):
    canonical_site: dict[int, CanonicalSite]
    xtal_form_site: dict[int, XtalFormSite]
    ligand_ids: list[LigandID]
    site_observation: dict[int, SiteObservation]

    @field_validator("canonical_site")
    @classmethod
    def check_canonical_site_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @field_validator("canonical_site")
    @classmethod
    def check_canonical_site_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @field_validator("xtal_form_site")
    @classmethod
    def check_xtal_form_site_ids(cls, v: dict[int, XtalFormSite]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @field_validator("xtal_form_site")
    @classmethod
    def check_xtal_form_site_ids_sequential(
        cls,
        v: dict[int, XtalFormSite],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @field_validator("site_observation")
    @classmethod
    def check_site_observation_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.ligand_id

    @field_validator("site_observation")
    @classmethod
    def check_site_observation_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums


def read_xmap(path: Path):
    m = gemmi.read_ccp4_map(str(path), setup=True)
    return m.grid


def transform_to_gemmi(transform: Transform):
    transform_gemmi = gemmi.Transform()
    transform_gemmi.vec.fromlist(transform.vec)
    transform_gemmi.mat.fromlist(transform.mat)

    return transform_gemmi


def gemmi_to_transform(transform):
    return Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist(), alignable_ids=[])


def get_box(neighbourhood: LigandNeighbourhood, xmap, transform):
    transform_gemmi = transform

    box = gemmi.FractionalBox()
    for atom in neighbourhood.atoms:
        transformed_pos = transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
        box.extend(
            xmap.unit_cell.fractionalize(gemmi.Position(transformed_pos.x, transformed_pos.y, transformed_pos.z))
        )

    for atom in neighbourhood.artefact_atoms:
        transformed_pos = transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
        box.extend(
            xmap.unit_cell.fractionalize(gemmi.Position(transformed_pos.x, transformed_pos.y, transformed_pos.z))
        )
    return box


def write_xmap(xmap, path: Path, neighbourhood: LigandNeighbourhood, transform):
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = xmap
    ccp4.setup(float("nan"))
    ccp4.update_ccp4_header()

    box = get_box(neighbourhood, xmap, transform)
    box_min = box.minimum
    box_max = box.maximum
    box_min_str = f"{round(box_min.x, 2)} {round(box_min.y, 2)} {round(box_min.z, 2)}"
    box_max_str = f"{round(box_max.x, 2)} {round(box_max.y, 2)} {round(box_max.z, 2)}"
    logger.debug(f"Box Extent is: min {box_min_str} : max {box_max_str}")
    ccp4.set_extent(box)
    ccp4.setup(float("nan"))
    ccp4.update_ccp4_header()

    ccp4.write_ccp4_map(str(path))


def read_graph(path: Path):
    g = nx.read_gml(
        str(path / constants.ALIGNABILITY_GRAPH_FILE_NAME),
        destringizer=lambda x: LigandID.from_string(x),
    )

    return g


def read_neighbourhoods(path: Path):
    neighbourhoods = LigandNeighbourhoods.parse_file(str(path / constants.NEIGHBOURHOODS_FILE_NAME))
    return neighbourhoods


def read_canonical_sites(path: Path):
    sites = CanonicalSites.parse_file(str(path / constants.SITES_FILE_NAME))

    return sites


def read_transforms(path: Path):
    transforms = Transforms.parse_file(str(path / constants.TRANSFORMS_FILE_NAME))
    return transforms


def read_structures(system_data: SystemData):
    structures = {}
    for dataset in system_data.datasets:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    return structures


def read_system_data(path: Path):
    return SystemData.parse_file(str(path / constants.DATA_JSON_PATH))


class SiteTransforms(BaseModel):
    canonical_site_transform_ids: list[tuple[int, int]]
    canonical_site_transforms: list[Transform]
    conformer_site_transform_ids: list[tuple[int, int, int]]
    conformer_site_transforms: list[Transform]

    def get_conformer_site_transform(self, site_id, subsite_id):
        for subsite_transform_id, subsite_transform in zip(
            self.conformer_site_transform_ids, self.conformer_site_transforms
        ):
            if subsite_transform_id[0] == site_id:
                if subsite_transform_id[2] == subsite_id:
                    return subsite_transform

        raise Exception()

    def get_canonical_site_transform(
        self,
        site_id,
    ):
        for site_transform_id, subsite_transform in zip(
            self.canonical_site_transform_ids, self.canonical_site_transforms
        ):
            if site_transform_id[1] == site_id:
                return subsite_transform

        raise Exception()


def save_site_transforms(site_transforms: SiteTransforms, path: Path):
    with open(path / constants.SITES_TRANSFORMS_FILE_NAME, "w") as f:
        f.write(site_transforms.json())


def read_site_transforms(path: Path):
    return SiteTransforms.parse_file(str(path / constants.SITES_TRANSFORMS_FILE_NAME))


def save_canonical_sites(sites: CanonicalSites, path: Path):
    with open(path / constants.SITES_FILE_NAME, "w") as f:
        f.write(sites.json())


def save_data(system_data: SystemData, output_dir: Path):
    with open(output_dir / constants.DATA_JSON_PATH, "w") as f:
        f.write(system_data.json())


class Options(BaseModel):
    source_dir: str
    output_dir: str
    datasources: list[str]
    datasource_types: list[str]
    panddas: list[str]
    assemblies_json: str
    xtalforms_json: str
    # dataset_xtalforms_json: str


class AssignedXtalForms(BaseModel):
    dataset_ids: list[DatasetID]
    xtalform_ids: list[int]

    def iter(self):
        for dataset_id, xtalform_id in zip(self.dataset_ids, self.xtalform_ids):
            yield dataset_id, xtalform_id

    def get_xtalform_id(self, item):
        for dataset_id, xtalform_id in self.iter():
            if dataset_id == item:
                return xtalform_id

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())


def read_assigned_xtalforms(path: Path):
    return AssignedXtalForms.parse_file(path / constants.ASSIGNED_XTALFORMS_FILE_NAME)


def save_assigned_xtalforms(path: Path, assigned_xtalforms: AssignedXtalForms):
    with open(path / constants.ASSIGNED_XTALFORMS_FILE_NAME, "w") as f:
        f.write(assigned_xtalforms.json())


def save_xtalforms(path: Path, xtalforms: XtalForms):
    with open(path / constants.XTALFORMS_FILE_NAME, "w") as f:
        f.write(xtalforms.json())


def read_xtalforms(path: Path):
    return XtalForms.parse_file(path / constants.XTALFORMS_FILE_NAME)


class LigandOutput(BaseModel):
    aligned_structures: dict[int, str]
    aligned_artefacts: dict[int, str]
    aligned_xmaps: dict[int, str]
    aligned_event_maps: dict[int, str]


class ChainOutput(BaseModel):
    aligned_ligands: dict[int, LigandOutput]

    def __getitem__(self, item):
        return self.aligned_ligands[item]


class DatasetOutput(BaseModel):
    aligned_chain_output: dict[str, ChainOutput]

    def __getitem__(self, item):
        return self.aligned_chain_output[item]


class Output(BaseModel):
    source_dir: str
    system_data: str
    xtalforms: str
    assigned_xtalforms: str
    neighbourhoods: str
    graph: str
    transforms: str
    sites: str
    site_transforms: str
    aligned_dir: str
    dataset_output: dict[str, DatasetOutput]

    @classmethod
    def read(cls, path: Path):
        return cls.parse_file(path)

    def save(self, path: Path):
        with open(path, "w") as f:
            f.write(self.json())

    def __getitem__(self, item):
        return self.dataset_output[item]


def save_output(
    output,
    path,
):
    with open(path / constants.OUTPUT_JSON_PATH, "w") as f:
        f.write(output.json())


def read_output(path):
    return Output.parse_file(path / constants.OUTPUT_JSON_PATH)


class LigandNeighbourhoodOutput:
    def __init__(
        self,
        aligned_structures: dict[str, str],
        aligned_artefacts: dict[str, str],
        aligned_xmaps: dict[str, str],
        aligned_diff_maps: dict[str, str],
        aligned_event_maps: dict[str, str],
        aligned_xmaps_crystallographic,
        aligned_diff_maps_crystallographic,
        aligned_event_maps_crystallographic,
    ):
        self.aligned_structures = aligned_structures
        self.aligned_artefacts: dict[str, str] = aligned_artefacts

        self.aligned_xmaps: dict[str, str] = aligned_xmaps
        self.aligned_diff_maps: dict[str, str] = aligned_diff_maps
        self.aligned_event_maps: dict[str, str] = aligned_event_maps

        self.aligned_xmaps_crystallographic: dict[str, str] = aligned_xmaps_crystallographic
        self.aligned_diff_maps_crystallographic: dict[str, str] = aligned_diff_maps_crystallographic
        self.aligned_event_maps_crystallographic: dict[str, str] = aligned_event_maps_crystallographic

    @staticmethod
    def from_dict(dic, source_dir):
        return LigandNeighbourhoodOutput(
            aligned_structures={k: Path(v) for k, v in dic["aligned_structures"].items()},
            aligned_artefacts={k: Path(v) for k, v in dic["aligned_artefacts"].items()},
            aligned_xmaps={k: Path(v) for k, v in dic["aligned_xmaps"].items()},
            aligned_diff_maps={k: Path(v) for k, v in dic["aligned_diff_maps"].items()},
            aligned_event_maps={k: Path(v) for k, v in dic["aligned_event_maps"].items()},
            aligned_xmaps_crystallographic={k: Path(v) for k, v in dic["aligned_xmaps_crystallographic"].items()},
            aligned_diff_maps_crystallographic={
                k: Path(v) for k, v in dic["aligned_diff_maps_crystallographic"].items()
            },
            aligned_event_maps_crystallographic={
                k: Path(v) for k, v in dic["aligned_event_maps_crystallographic"].items()
            },
        )

    def to_dict(self):
        dic = {
            "aligned_structures": {
                canonical_site_id: str(path) for canonical_site_id, path in self.aligned_structures.items()
            },
            "aligned_artefacts": {
                canonical_site_id: str(path) for canonical_site_id, path in self.aligned_artefacts.items()
            },
            "aligned_xmaps": {canonical_site_id: str(path) for canonical_site_id, path in self.aligned_xmaps.items()},
            "aligned_diff_maps": {
                canonical_site_id: str(path) for canonical_site_id, path in self.aligned_diff_maps.items()
            },
            "aligned_event_maps": {
                canonical_site_id: str(path) for canonical_site_id, path in self.aligned_event_maps.items()
            },
            "aligned_xmaps_crystallographic": {
                canonical_site_id: str(path) for canonical_site_id, path in self.aligned_xmaps_crystallographic.items()
            },
            "aligned_diff_maps_crystallographic": {
                canonical_site_id: str(path)
                for canonical_site_id, path in self.aligned_diff_maps_crystallographic.items()
            },
            "aligned_event_maps_crystallographic": {
                canonical_site_id: str(path)
                for canonical_site_id, path in self.aligned_event_maps_crystallographic.items()
            },
        }
        return dic


def symlink(old_path, new_path):
    os.symlink(old_path.resolve(), new_path)


class FSModel:
    def __init__(
        self,
        source_dir,
        fs_model,
        # assemblies,
        xtalforms,
        dataset_assignments,
        ligand_neighbourhoods,
        alignability_graph,
        connected_components,
        ligand_neighbourhood_transforms,
        conformer_sites,
        conformer_site_transforms,
        canonical_sites,
        xtalform_sites,
        reference_structure_transforms,
        alignments,
        reference_alignments,
        hierarchy,
        biochain_priorities,
        assembly_landmarks,
        assembly_transforms,
        chain_to_assembly,
    ):
        self.source_dir = source_dir
        self.fs_model = fs_model
        self.xtalforms = xtalforms
        self.dataset_assignments = dataset_assignments
        self.ligand_neighbourhoods = ligand_neighbourhoods
        self.alignability_graph = alignability_graph
        self.connected_components = connected_components
        self.ligand_neighbourhood_transforms = ligand_neighbourhood_transforms
        self.conformer_sites = conformer_sites
        self.conformer_site_transforms = conformer_site_transforms
        self.canonical_sites = canonical_sites
        self.xtalform_sites = xtalform_sites
        self.reference_structure_transforms = reference_structure_transforms
        self.alignments = alignments
        self.reference_alignments = reference_alignments

        self.hierarchy = hierarchy
        self.biochain_priorities = biochain_priorities
        self.assembly_landmarks = assembly_landmarks
        self.assembly_transforms = assembly_transforms
        self.chain_to_assembly = chain_to_assembly

    def symlink_old_data(self):
        for dtag, dataset_alignments in self.alignments.items():
            if not (self.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                os.mkdir(self.source_dir / constants.ALIGNED_FILES_DIR / dtag)
            for chain, chain_alignments in dataset_alignments.items():
                for residue, ligand_neighbourhood_alignments in chain_alignments.items():
                    for canonical_site_id in ligand_neighbourhood_alignments.aligned_structures:
                        if canonical_site_id in ligand_neighbourhood_alignments.aligned_structures:
                            old_path = Path(ligand_neighbourhood_alignments.aligned_structures[canonical_site_id])
                            new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                            if not new_path.exists():
                                symlink(old_path, new_path)
                        if canonical_site_id in ligand_neighbourhood_alignments.aligned_artefacts:
                            old_path = Path(ligand_neighbourhood_alignments.aligned_artefacts[canonical_site_id])
                            new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                            if not new_path.exists():
                                symlink(old_path, new_path)
                        if canonical_site_id in ligand_neighbourhood_alignments.aligned_xmaps:
                            old_path = Path(ligand_neighbourhood_alignments.aligned_xmaps[canonical_site_id])
                            new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                            if not new_path.exists():
                                symlink(old_path, new_path)
                        if canonical_site_id in ligand_neighbourhood_alignments.aligned_event_maps:
                            old_path = Path(ligand_neighbourhood_alignments.aligned_event_maps[canonical_site_id])
                            new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                            if not new_path.exists():
                                symlink(old_path, new_path)

        # Symlink old alignments
        for dtag, dtag_alignment_info in self.reference_alignments.items():
            if not (self.source_dir / constants.ALIGNED_FILES_DIR / dtag).exists():
                os.mkdir(self.source_dir / constants.ALIGNED_FILES_DIR / dtag)
            for canonical_site_id, canonical_site_alignment_info in dtag_alignment_info.items():
                old_path = Path(canonical_site_alignment_info["aligned_structures"])
                new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                if not new_path.exists():
                    symlink(old_path, new_path)

                old_path = Path(canonical_site_alignment_info["aligned_artefacts"])
                new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                if not new_path.exists():
                    symlink(old_path, new_path)

                old_path = Path(canonical_site_alignment_info["aligned_xmaps"])
                new_path = self.source_dir / constants.ALIGNED_FILES_DIR / dtag / old_path.name
                if not new_path.exists():
                    symlink(old_path, new_path)

    @staticmethod
    def from_dir(
        source_dir: str,
        #        output_dir: str,
    ):
        source_dir = Path(source_dir)
        # output_dir = Path(output_dir)

        fs_model = source_dir / constants.FS_MODEL_YAML_FILE_NAME
        if fs_model.exists():
            with open(fs_model, "r") as f:
                dic = yaml.safe_load(f)
            if dic is not None:
                return FSModel.from_dict(dic)

        xtalforms = source_dir / constants.XTALFORMS_YAML_FILE_NAME
        dataset_assignments = source_dir / constants.ASSIGNED_XTALFORMS_YAML_FILE_NAME
        ligand_neighbourhoods = source_dir / constants.NEIGHBOURHOODS_YAML_FILE_NAME
        alignability_graph = source_dir / constants.ALIGNABILITY_GRAPH_FILE_NAME
        connected_components = source_dir / constants.CONNECTED_COMPONENTS_YAML_NAME
        ligand_neighbourhood_transforms = source_dir / constants.TRANSFORMS_YAML_FILE_NAME
        conformer_sites = source_dir / constants.CONFORMER_SITE_YAML_FILE
        conformer_site_transforms = source_dir / constants.CONFORMER_SITES_TRANSFORMS_YAML_FILE_NAME
        canonical_sites = source_dir / constants.CANONICAL_SITE_YAML_FILE
        xtalform_sites = source_dir / constants.XTALFORM_SITE_YAML_FILE
        reference_structure_transforms = source_dir / constants.REFERENCE_STRUCTURE_TRANSFORMS_YAML
        alignments = {}
        reference_alignments = {}

        hierarchy = source_dir / constants.HIERARCHY_YAML
        biochain_priorities = source_dir / constants.BIOCHAIN_PRIORITIES_YAML
        assembly_landmarks = source_dir / constants.ASSEMBLY_LANDMARKS_YAML
        assembly_transforms = source_dir / constants.ASSEMBLY_TRANSFORMS_YAML
        chain_to_assembly = source_dir / constants.CHAIN_TO_ASSEMBLY_YAML

        return FSModel(
            source_dir,
            fs_model,
            xtalforms,
            dataset_assignments,
            ligand_neighbourhoods,
            alignability_graph,
            connected_components,
            ligand_neighbourhood_transforms,
            conformer_sites,
            conformer_site_transforms,
            canonical_sites,
            xtalform_sites,
            reference_structure_transforms,
            alignments,
            reference_alignments,
            hierarchy,
            biochain_priorities,
            assembly_landmarks,
            assembly_transforms,
            chain_to_assembly,
        )

    @staticmethod
    def from_dict(dic):
        source_dir = Path(dic["source_dir"])
        alignments = {}
        for dtag, dataset_alignments in dic["alignments"].items():
            alignments[dtag] = {}
            for chain, chain_alignments in dataset_alignments.items():
                alignments[dtag][chain] = {}
                for residue, residue_alignments in chain_alignments.items():
                    alignments[dtag][chain][residue] = {}
                    for altloc, altloc_alignments in residue_alignments.items():
                        altloc = string_to_altloc(altloc)
                        alignments[dtag][chain][residue][altloc] = {}
                        for version, ligand_neighbourhood_alignments in altloc_alignments.items():
                            alignments[dtag][chain][residue][altloc][version] = LigandNeighbourhoodOutput.from_dict(
                                ligand_neighbourhood_alignments, source_dir
                            )

        reference_alignments = {}
        for dtag, canonical_site_alignments in dic["reference_alignments"].items():
            reference_alignments[dtag] = {}
            for canonical_site_id, canonical_site_alignment_info in canonical_site_alignments.items():
                reference_alignments[dtag][canonical_site_id] = {
                    "aligned_structures": Path(canonical_site_alignment_info["aligned_structures"]),
                    "aligned_artefacts": Path(canonical_site_alignment_info["aligned_artefacts"]),
                    "aligned_xmaps": Path(canonical_site_alignment_info["aligned_xmaps"]),
                }

        return FSModel(
            source_dir=Path(dic["source_dir"]),
            fs_model=Path(dic["fs_model"]),
            xtalforms=Path(dic["crystalforms"]),
            dataset_assignments=Path(dic["dataset_assignments"]),
            ligand_neighbourhoods=Path(dic["ligand_neighbourhoods"]),
            alignability_graph=Path(dic["alignability_graph"]),
            connected_components=Path(dic["connected_components"]),
            ligand_neighbourhood_transforms=Path(dic["ligand_neighbourhood_transforms"]),
            conformer_sites=Path(dic["conformer_sites"]),
            conformer_site_transforms=Path(dic["conformer_site_transforms"]),
            canonical_sites=Path(dic["canonical_sites"]),
            xtalform_sites=Path(dic["xtalform_sites"]),
            reference_structure_transforms=Path(dic["reference_structure_transforms"]),
            alignments=alignments,
            reference_alignments=reference_alignments,
            hierarchy=Path(dic["hierarchy"]),
            biochain_priorities=Path(dic["biochain_priorities"]),
            assembly_landmarks=Path(dic["assembly_landmarks"]),
            assembly_transforms=Path(dic["assembly_transforms"]),
            chain_to_assembly=Path(dic["chain_to_assembly"]),
        )

    def to_dict(
        self,
    ):
        dic = {}
        alignments = {}
        for dtag, dataset_alignments in self.alignments.items():
            alignments[dtag] = {}
            for chain, chain_alignments in dataset_alignments.items():
                alignments[dtag][chain] = {}
                for residue, residue_alignments in chain_alignments.items():
                    alignments[dtag][chain][residue] = {}
                    for altloc, altloc_alignments in residue_alignments.items():
                        altloc = altloc_to_string(altloc)
                        alignments[dtag][chain][residue][altloc] = {}
                        for version, ligand_neighbourhood_alignments in altloc_alignments.items():
                            alignments[dtag][chain][residue][altloc][version] = LigandNeighbourhoodOutput.to_dict(
                                ligand_neighbourhood_alignments
                            )

        reference_alignments = {}
        for dtag, dtag_alignment_info in self.reference_alignments.items():
            reference_alignments[dtag] = {}
            for canonical_site_id, canonical_site_alignment_info in dtag_alignment_info.items():
                reference_alignments[dtag][canonical_site_id] = {
                    "aligned_structures": str(canonical_site_alignment_info["aligned_structures"]),
                    "aligned_artefacts": str(canonical_site_alignment_info["aligned_artefacts"]),
                    "aligned_xmaps": str(canonical_site_alignment_info["aligned_xmaps"]),
                }

        return {
            "source_dir": str(self.source_dir),
            "fs_model": str(self.fs_model),
            "crystalforms": str(self.xtalforms),
            "dataset_assignments": str(self.dataset_assignments),
            "ligand_neighbourhoods": str(self.ligand_neighbourhoods),
            "alignability_graph": str(self.alignability_graph),
            "connected_components": str(self.connected_components),
            "ligand_neighbourhood_transforms": str(self.ligand_neighbourhood_transforms),
            "conformer_sites": str(self.conformer_sites),
            "conformer_site_transforms": str(self.conformer_site_transforms),
            "canonical_sites": str(self.canonical_sites),
            "xtalform_sites": str(self.xtalform_sites),
            "reference_structure_transforms": str(self.reference_structure_transforms),
            "alignments": alignments,
            "reference_alignments": reference_alignments,
            'hierarchy': str(self.hierarchy),
            'biochain_priorities': str(self.biochain_priorities),
            'assembly_landmarks': str(self.assembly_landmarks),
            'assembly_transforms': str(self.assembly_transforms),
            'chain_to_assembly': str(self.chain_to_assembly),
        }


class Datasource:
    def __init__(self, path: str, datasource_type: str):
        self.path = path
        self.datasource_type = datasource_type


class PanDDA:
    def __init__(
        self,
        path: str,
    ):
        self.path = Path(path)
        self.event_table_path = self.path / constants.PANDDA_ANALYSES_DIR / constants.PANDDA_EVENTS_INSPECT_TABLE_PATH


class LigandBindingEvent:
    def __init__(
        self,
        id,
        dtag,
        chain,
        residue,
        altloc,
        xmap,
    ):
        self.id: str = id
        self.dtag: str = dtag
        self.chain: str = chain
        self.residue: str = residue
        self.altloc:str = altloc
        self.xmap: str = xmap




# class SourceDataModel:
#     def __init__(
#         self,
#         fs_model: FSModel,
#         datasources: list[Datasource],
#         panddas: list[PanDDA],
#     ):
#         self.fs_model = fs_model
#         self.datasources = datasources
#         self.panddas = panddas

#     @staticmethod
#     def from_fs_model(
#         fs_model: FSModel,
#         datasources,
#         datasource_types,
#         panddas,
#     ):
#         _datasources = []
#         for datasource_path, datasource_type in zip(datasources, datasource_types):
#             datasource = Datasource(datasource_path, datasource_type)
#             _datasources.append(datasource)

#         _panddas = []
#         for pandda_dir in panddas:
#             pandda = PanDDA(pandda_dir)
#             _panddas.append(pandda)

#         return SourceDataModel(fs_model, _datasources, _panddas)

#     def get_datasets(self):
#         datasets = {}
#         reference_datasets = {}
#         new_datasets = {}

#         # Get the pandda tables
#         pandda_event_tables = {pandda.path: pd.read_csv(pandda.event_table_path) for pandda in self.panddas}

#         # Get all the datasets attested in the data sources
#         for datasource in self.datasources:
#             logger.info(f"Parsing datasource: {datasource.path}")
#             if datasource.datasource_type == "model_building":
#                 for model_dir in Path(datasource.path).glob("*"):
#                     dtag = model_dir.name
#                     if dtag in datasets:
#                         st = f"Dataset ID {dtag} already found! Using new!"
#                         logger.warning(st)
#                         continue

#                     pdb = model_dir / constants.MODEL_DIR_PDB
#                     xmap = model_dir / constants.MODEL_DIR_XMAP
#                     mtz = model_dir / constants.MODEL_DIR_MTZ
#                     if not pdb.exists():
#                         continue

#                     ligand_binding_events = _get_ligand_binding_events_from_panddas(
#                         pandda_event_tables,
#                         pdb,
#                         dtag,
#                     )
#                     if len(ligand_binding_events) == 0:
#                         logger.warning(f"Dataset {dtag} has no ligand binding events!")
#                         continue
#                     dataset = Dataset(
#                         dtag=dtag,
#                         pdb=str(pdb),
#                         xmap=str(xmap),
#                         mtz=str(mtz),
#                         ligand_binding_events=ligand_binding_events,
#                     )
#                     datasets[dtag] = dataset
#                     logger.debug(f"Added dataset: {dtag}")

#             elif datasource.datasource_type == "manual":
#                 for model_dir in Path(datasource.path).glob("*"):
#                     dtag = model_dir.name

#                     if dtag in datasets:
#                         st = f"Dataset ID {dtag} already found! Using new!"
#                         logger.warning(st)
#                     try:
#                         pdb = next(model_dir.glob("*.pdb"))
#                     except Exception:
#                         raise Exception(f"Could not find pdb in dir: {model_dir}")
#                     try:
#                         xmap = next(model_dir.glob("*.ccp4"))
#                     except Exception as e:
#                         print(e)
#                         xmap = None
#                         logger.warning("No xmap!")
#                     try:
#                         mtz = next(model_dir.glob("*.mtz"))
#                     except Exception as e:
#                         print(e)
#                         mtz = None
#                         logger.warning("No mtz!")

#                     ligand_binding_events = _get_ligand_binding_events_from_structure(pdb, xmap, dtag)
#                     if len(ligand_binding_events) == 0:
#                         logger.warning(f"Dataset {dtag} has no ligand binding events!")
#                         continue
#                     dataset = Dataset(
#                         dtag=dtag,
#                         pdb=str(pdb),
#                         xmap=str(xmap),
#                         mtz=str(mtz),
#                         ligand_binding_events=ligand_binding_events,
#                     )
#                     datasets[dtag] = dataset
#                     reference_datasets[dtag] = dataset
#                     logger.debug(f"Added dataset: {dtag}")
#             else:
#                 raise Exception(f"Source type {datasource.datasource_type} unknown!")

#         # Determine which of these are new using the fs_model output
#         for dtag, dataset in datasets.items():
#             if dtag not in self.fs_model.alignments:
#                 new_datasets[dtag] = dataset

#         return datasets, reference_datasets, new_datasets

#     @staticmethod
#     def from_dict():
#         ...

#     def to_dict(self, path: Path):
#         ...


class Generator:
    def __init__(self, biomol: str, chain: str, triplet: str):
        self.biomol: str = biomol
        self.chain: str = chain
        self.triplet: str = triplet


class Assembly:
    def __init__(self, reference: str, generators: list[Generator]):
        self.reference = reference
        self.generators = generators

    @staticmethod
    def from_dict(dic):
        reference = dic["reference"]
        biomol = dic["biomol"]
        chains = dic["chains"]

        # Split biomol on commas and strip whitespace
        biomol_matches = re.findall("([A-Z]+)", biomol)

        # Split chains on commas that do not follow a number, x,y or z and strip whitespace
        chain_matches = re.findall("(([A-Z]+)([(]+[^()]+[)]+)*)", chains)
        # print("biomol_matches:", biomol_matches)
        # print("chain_matches:", chain_matches)
        # Make generators
        generators = []
        for biomol_match, chain_match in zip(biomol_matches, chain_matches):
            if len(chain_match[2]) == 0:
                xyz = "x,y,z"
            else:
                xyz = chain_match[2][1:-1]
            generators.append(Generator(biomol_match, chain_match[1], xyz))

        return Assembly(reference, generators)


class XtalFormAssembly:
    def __init__(self, assembly: str, chains: list[str], transforms: list[str]):
        self.assembly = assembly
        self.chains = chains
        self.transforms = transforms


class XtalForm:
    def __init__(self, reference: str, assemblies: dict[str, XtalFormAssembly]):
        self.reference = reference
        self.assemblies = assemblies

    @staticmethod
    def from_dict(dic):
        reference = dic["reference"]
        assemblies = dic["assemblies"]

        _assemblies = {}
        for xtalform_assembly_id, xtalform_assembly_info in assemblies.items():
            assemblies[xtalform_assembly_id] = {}
            assembly = xtalform_assembly_info["assembly"]
            chains = xtalform_assembly_info["chains"]
            chains_matches = re.findall("(([A-Z]+)([(]+[^()]+[)]+)*)", chains)
            _chains = []
            _transforms = []
            for chain_match in chains_matches:
                _chains.append(chain_match[1])

                if len(chain_match[2]) == 0:
                    xyz = "x,y,z"
                else:
                    xyz = chain_match[2][1:-1]
                _transforms.append(xyz)

            _assemblies[xtalform_assembly_id] = XtalFormAssembly(assembly, _chains, _transforms)
        return XtalForm(reference, _assemblies)


def altloc_to_string(altloc):
    if altloc == '\0':
        return '0'
    else:
        return altloc

def string_to_altloc(altloc):
    if str(altloc) == '0':
        return '\0'
    else:
        return altloc 


class Transform:
    def __init__(self, vec, mat, alignable_ids):
        self.vec: list[float] = vec
        self.mat: list[list[float]] = mat
        self.alignable_ids = alignable_ids


    @staticmethod
    def from_dict(dic):
        return Transform(dic["vec"], dic["mat"], dic['alignable_ids'])

    def to_dict(self):
        return {"vec": self.vec, "mat": self.mat, "alignable_ids": self.alignable_ids}


def transform_to_string(transform):
    return str([x for x in transform.vec] + [x for y in transform.mat for x in y])


class Atom:
    def __init__(
        self,
        element: str,
        x: float,
        y: float,
        z: float,
        image: Transform,
    ):
        self.element: str = element
        self.x: float = x
        self.y: float = y
        self.z: float = z
        self.image: Transform = image

    @staticmethod
    def from_dict(dic):
        return Atom(dic["element"], dic["x"], dic["y"], dic["z"], Transform.from_dict(dic["image"]))

    def to_dict(self):
        dic = {}
        return {"element": self.element, "x": self.x, "y": self.y, "z": self.z, "image": self.image.to_dict()}


class Neighbourhood:
    def __init__(self, atoms: dict[tuple[str, str, str], Atom], artefact_atoms: dict[tuple[str, str, str], Atom]):
        self.atoms = atoms
        self.artefact_atoms = artefact_atoms

    @staticmethod
    def from_dict(dic):
        atoms = {}
        artefact_atoms = {}

        _atoms = dic["atoms"]
        _artefact_atoms = dic["artefact_atoms"]
        for atom_id, atom_info in _atoms.items():
            chain, residue, atom = atom_id.split("/")
            atoms[(chain, residue, atom)] = Atom.from_dict(atom_info)
            ...

        for atom_id, atom_info in _artefact_atoms.items():
            chain, residue, atom = atom_id.split("/")
            artefact_atoms[(chain, residue, atom)] = Atom.from_dict(atom_info)

        return Neighbourhood(atoms, artefact_atoms)

    def to_dict(self):
        dic = {}
        dic["atoms"] = {}
        for atom_id, atom in self.atoms.items():
            dic["atoms"]["/".join(atom_id)] = atom.to_dict()
        dic["artefact_atoms"] = {}
        for atom_id, atom in self.artefact_atoms.items():
            dic["artefact_atoms"]["/".join(atom_id)] = atom.to_dict()

        return dic


class AlignabilityGraph:
    ...


class ConformerSite:
    def __init__(
        self,
        residues: list[tuple[str, str]],
        residues_aligned: list[tuple[str, str]],
        members: list[tuple[str, str, str, str]],
        reference_ligand_id: tuple[str, str, str, str],
    ):
        self.residues: list[tuple[str, str]] = residues
        self.residues_aligned = residues_aligned
        self.members: list[tuple[str, str, str, str]] = members
        self.reference_ligand_id: tuple[str, str, str, str] = reference_ligand_id

    @staticmethod
    def from_dict(dic):
        residues = []
        for res in dic["residues"]:
            chain, residue, name = res.split("/")
            residues.append((chain, residue, name))
        residues_aligned = []
        for res in dic["residues_aligned"]:
            chain, residue, name = res.split("/")
            residues.append((chain, residue, name))

        members = []
        for member in dic["members"]:
            dtag, chain, residue, altloc, version = [string_to_altloc(x) for x in member.split("/")]
            members.append((dtag, chain, residue, string_to_altloc(altloc), version))
        ref_dtag, ref_chain, ref_residue, ref_altloc, version = [string_to_altloc(x) for x in dic["reference_ligand_id"].split("/")]
        return ConformerSite(residues, residues_aligned, members, (ref_dtag, ref_chain, ref_residue, ref_altloc, version))

    def to_dict(
        self,
    ):
        return {
            "residues": [x for x in sorted(["/".join(resid) for resid in self.residues])],
            "residues_aligned": [x for x in sorted(["/".join(resid) for resid in self.residues_aligned])],
            "members": [x for x in sorted(["/".join([altloc_to_string(x) for x in lid]) for lid in self.members])],
            "reference_ligand_id": "/".join([altloc_to_string(x) for x in self.reference_ligand_id]),
        }


class CanonicalSite:
    def __init__(
        self,
        conformer_site_ids: list[str],
        residues: list[tuple[str, str]],
        reference_conformer_site_id: str,
        global_reference_dtag: str,
        centroid_res: tuple[str, str, str, str],
    ):
        self.conformer_site_ids: list[str] = conformer_site_ids
        self.residues: list[tuple[str, str]] = residues
        self.reference_conformer_site_id: str = reference_conformer_site_id
        self.global_reference_dtag: str = global_reference_dtag
        self.centroid_res = centroid_res

    @staticmethod
    def from_dict(dic):
        residues = []
        for res in dic["residues"]:
            chain, residue, name = res.split("/")
            residues.append((chain, residue, name))

        return CanonicalSite(
            dic["conformer_site_ids"],
            residues,
            dic["reference_conformer_site_id"],
            dic["global_reference_dtag"],
            dic["centroid_res"].split("/"),
        )

    def to_dict(self):
        return {
            "conformer_site_ids": [altloc_to_string(x) for x in self.conformer_site_ids],
            "residues": [x for x in sorted(["/".join(res) for res in self.residues])],
            "reference_conformer_site_id": self.reference_conformer_site_id,
            "global_reference_dtag": self.global_reference_dtag,
            "centroid_res": "/".join(self.centroid_res),
        }


class XtalFormSite:
    def __init__(
        self,
        xtalform_id: str,
        crystallographic_chain: str,
        canonical_site_id: str,
        members: list[tuple[str, str, str]],
    ):
        self.xtalform_id: str = xtalform_id
        self.crystallographic_chain: str = crystallographic_chain
        self.canonical_site_id: str = canonical_site_id
        self.members: list[tuple[str, str, str]] = members

    @staticmethod
    def from_dict(dic):
        members = []
        for member in dic["members"]:
            dtag, chain, residue, altloc, version = [string_to_altloc(x) for x in member.split("/")]
            members.append((dtag, chain, residue, string_to_altloc(altloc), version))
        return XtalFormSite(dic["xtalform_id"], dic["crystallographic_chain"], dic["canonical_site_id"], members)

    def to_dict(self):
        dic = {}
        dic["members"] = []
        for member in self.members:
            dic["members"].append("/".join([altloc_to_string(x) for x in member]))

        dic["xtalform_id"] = self.xtalform_id
        dic["crystallographic_chain"] = self.crystallographic_chain
        dic["canonical_site_id"] = self.canonical_site_id
        return dic

    def __rich_repr__(self):
        yield self.to_dict()



class Dataset:
    def __init__(
        self,
        dtag,
        pdb,
        xmap,
        mtz,
        ligand_binding_events: dict[LigandNeighbourhoodID, LigandBindingEvent],
    ):
        self.dtag = dtag
        self.pdb = pdb
        self.xmap = xmap
        self.mtz = mtz
        self.ligand_binding_events = ligand_binding_events

