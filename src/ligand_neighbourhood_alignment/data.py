from pathlib import Path
from typing import Generator

import gemmi
import networkx as nx
from loguru import logger
from pydantic import BaseModel, field_validator

from ligand_neighbourhood_alignment import constants


Structure = gemmi.Structure


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


# class CanonicalSite(BaseModel):
#     id: int
#     name: str
#     refpdb: str
#     atoms: dict[AtomID, Atom]
#     literatureref: str
#     members: list[LigandID]


# class XtalForm(BaseModel):
#     id: int
#     space_group: int
#     unit_cell: tuple[float, float, float, float, float, float]
#     members: list[DatasetID]
#     transforms: list[Transform]


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
    # xtalform_ids: list[int]
    # xtalforms: list[XtalForm]
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


# class XtalFormSite(BaseModel):
#     id: int
#     canon_site_id: int
#     xtal_form_id: int
#     code: str
#     refpdb: str
#     atoms: dict[AtomID, Atom]
#     artefact_atoms: dict[AtomID, Atom]
#     members: list[LigandID]


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
    # data_source_type: str
    datasource_type: str
    # dataset_ids: list[DatasetID]
    # datasets: list[Dataset]


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


# class Site(BaseModel):
#     members: list[LigandID]


class Block(BaseModel):
    xi: int
    yi: int
    zi: int
    xmi: int
    ymi: int
    zmi: int
    dx: int
    dy: int
    dz: int
    transform: Transform


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
    # subsites: list[ConformerSite]
    # reference_subsite: int
    conformer_sites: dict[int, ConformerSite]

    # def iter(self) -> Generator[tuple[int, ConformerSite], None, None]:
    #     for subsite in self.subsites:
    #         yield subsite.id, subsite

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
    return Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist())


def get_box(neighbourhood: LigandNeighbourhood, xmap, transform):
    # transform_gemmi = transform_to_gemmi(transform)
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
    # aligned_structures: dict[str, ChainOutput]
    # aligned_xmaps: dict[str, ChainOutput]

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
