import yaml
import gemmi
import networkx as nx

from ligand_neighbourhood_alignment import dt

def _get_structures(datasets):
    structures = {}
    for dtag, dataset in datasets.items():
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dtag] = structure

    return structures


def _save_neighbourhoods(
    fs_model: dt.FSModel,
    ligand_neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
):
    with open(fs_model.ligand_neighbourhoods, "w") as f:
        dic = {}
        for ligand_id, neighbourhood in ligand_neighbourhoods.items():
            dic["/".join([dt.altloc_to_string(x) for x in ligand_id])] = neighbourhood.to_dict()
        yaml.safe_dump(dic, f)


def _save_ligand_neighbourhood_transforms(fs_model, ligand_neighbourhood_transforms):
    with open(fs_model.ligand_neighbourhood_transforms, "w") as f:
        dic = {}
        for (to_ligand_id, from_ligand_id), transform in ligand_neighbourhood_transforms.items():
            key = "~".join(["/".join(to_ligand_id), "/".join(from_ligand_id)])
            dic[key] = transform.to_dict()
        yaml.safe_dump(dic, f)


def _save_assignments(fs_model: dt.FSModel, dataset_assignments: dict[str, str]):
    with open(fs_model.dataset_assignments, "w") as f:
        yaml.safe_dump(dataset_assignments, f)

def _save_graph(fs_model, alignability_graph):
    graph_for_output = nx.relabel_nodes(alignability_graph, {x: "/".join(x) for x in alignability_graph})
    nx.write_gml(
        graph_for_output,
        str(fs_model.alignability_graph),
        # stringizer=lambda x: "/".join(x),
    )

def _save_connected_components(fs_model, connected_components):
    with open(fs_model.connected_components, "w") as f:
        dic = {}
        for connected_component_reference, connected_component in connected_components.items():
            dic["+".join(connected_component_reference)] = ["+".join(member) for member in connected_component]
        yaml.safe_dump(dic, f, sort_keys=False)


def _save_conformer_sites(fs_model: dt.FSModel, conformer_sites: dict[str, dt.ConformerSite]):
    with open(fs_model.conformer_sites, "w") as f:
        dic = {}
        for conformer_site_id, conformer_site in conformer_sites.items():
            dic[dt.altloc_to_string(conformer_site_id)] = conformer_site.to_dict()
        yaml.safe_dump(dic, f)


def _save_canonical_sites(fs_model, canonical_sites: dict[str, dt.CanonicalSite]):
    with open(fs_model.canonical_sites, "w") as f:
        dic = {}
        for canonical_site_id, canonical_site in canonical_sites.items():
            dic[dt.altloc_to_string(canonical_site_id)] = canonical_site.to_dict()
        yaml.safe_dump(dic, f)

def _save_xtalform_sites(fs_model, xtalform_sites: dict[str, dt.XtalFormSite]):
    with open(fs_model.xtalform_sites, "w") as f:
        dic = {}
        for xtalform_site_id, xtalform_site in xtalform_sites.items():
            dic[dt.altloc_to_string(xtalform_site_id)] = xtalform_site.to_dict()
        yaml.safe_dump(dic, f)


def _save_conformer_site_transforms(
    fs_model: dt.FSModel, conformer_site_transforms: dict[tuple[str, str], dt.Transform]
):
    with open(fs_model.conformer_site_transforms, "w") as f:
        dic = {}
        for conformer_site_transform_id, conformer_site_transform in conformer_site_transforms.items():
            dic["~".join(conformer_site_transform_id)] = conformer_site_transform.to_dict()
        yaml.safe_dump(dic, f)
    ...


def _save_canonical_site_transforms(fs_model: dt.FSModel, canonical_site_transforms: dict[str, dt.Transform]):
    with open(fs_model.canonical_site_transforms, "w") as f:
        dic = {}
        for canonical_site_transform_id, canonical_site_transform in canonical_site_transforms.items():
            dic[canonical_site_transform_id] = canonical_site_transform.to_dict()
        yaml.safe_dump(dic, f)

def _save_fs_model(fs_model: dt.FSModel):
    with open(fs_model.fs_model, "w") as f:
        dic = fs_model.to_dict()

        yaml.safe_dump(dic, f)


def _save_reference_structure_transforms(
    fs_model: dt.FSModel,
    reference_structure_transforms: dict[tuple[str, str], dt.Transform],
):
    dic = {}
    for reference_structure_transform_id, reference_structure_transform in reference_structure_transforms.items():
        dic["~".join(reference_structure_transform_id)] = reference_structure_transform.to_dict()
    with open(fs_model.reference_structure_transforms, "w") as f:
        yaml.safe_dump(dic, f)

def save_yaml(path, obj, obj_to_dict):
    with open(path, "w") as f:
        dic = obj_to_dict(obj)
        yaml.safe_dump(dic, f, sort_keys=False)


def load_yaml(path, dict_to_obj):
    with open(path, "r") as f:
        dic = yaml.safe_load(f)

    return dict_to_obj(dic)


def _load_assemblies(assemblies_file, new_assemblies_yaml):
    assemblies = {}

    if assemblies_file.exists():
        with open(assemblies_file, "r") as f:
            dic = yaml.safe_load(f)

        for assembly_id, assembly_info in dic["assemblies"].items():
            assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    # Load new info and update
    if new_assemblies_yaml.exists():
        with open(new_assemblies_yaml, "r") as f:
            new_assemblies_dict = yaml.safe_load(f)["assemblies"]
    else:
        new_assemblies_dict = {}

    for assembly_id, assembly_info in new_assemblies_dict.items():
        if assembly_id in assemblies:
            continue
        assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return assemblies


def _load_xtalforms(xtalforms_file, new_xtalforms_yaml):
    xtalforms = {}

    if xtalforms_file.exists():
        with open(xtalforms_file, "r") as f:
            dic = yaml.safe_load(f)["crystalforms"]

        for xtalform_id, xtalform_info in dic.items():
            xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    # Load new info and update
    if new_xtalforms_yaml.exists():
        with open(new_xtalforms_yaml, "r") as f:
            new_xtalforms_dict = yaml.safe_load(f)["crystalforms"]
    else:
        new_xtalforms_dict = {}

    for xtalform_id, xtalform_info in new_xtalforms_dict.items():
        if xtalform_id in xtalforms:
            continue
        xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return xtalforms


def _load_xtalforms_and_assemblies(xtalforms_file, new_xtalforms_yaml):
    assemblies = {}
    xtalforms = {}

    if xtalforms_file.exists():
        with open(xtalforms_file, "r") as f:
            dic = yaml.safe_load(f)

        for xtalform_id, xtalform_info in dic.items():
            xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    # Load new info and update
    if new_xtalforms_yaml.exists():
        with open(new_xtalforms_yaml, "r") as f:
            new_xtalforms_dict = yaml.safe_load(f)
    else:
        new_xtalforms_dict = {}

    for xtalform_id, xtalform_info in new_xtalforms_dict.items():
        if xtalform_id in xtalforms:
            continue
        xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return xtalforms


def _load_dataset_assignments(dataset_assignments_yaml, fail_if_not_found=False):
    if fail_if_not_found and not dataset_assignments_yaml.is_file():
        raise ValueError("File " + str(dataset_assignments_yaml) + " does not exist")
    dataset_assignments = {}
    if dataset_assignments_yaml.exists():
        with open(dataset_assignments_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for dtag, assignment in dic.items():
            dataset_assignments[dtag] = assignment

    return dataset_assignments


def _load_ligand_neighbourhoods(ligand_neighbourhoods_yaml, fail_if_not_found=False):
    if fail_if_not_found and not ligand_neighbourhoods_yaml.is_file():
        raise ValueError("File " + str(ligand_neighbourhoods_yaml) + " does not exist")
    ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = {}

    if ligand_neighbourhoods_yaml.exists():
        with open(ligand_neighbourhoods_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic:
            for ligand_id, neighbourhood_info in dic.items():
                dtag, chain, residue, altloc, version = [dt.string_to_altloc(x) for x in ligand_id.split("/")]
                neighbourhood = dt.Neighbourhood.from_dict(neighbourhood_info)
                ligand_neighbourhoods[(dtag, chain, residue, altloc, version)] = neighbourhood

    return ligand_neighbourhoods


def _load_alignability_graph(alignability_graph, fail_if_not_found=False):
    if fail_if_not_found and not alignability_graph.is_file():
        raise ValueError("File " + str(alignability_graph) + " does not exist")
    if alignability_graph.exists():
        g_initial = nx.read_gml(
            str(alignability_graph),
        )

        g = nx.relabel_nodes(g_initial, {x: tuple([dt.string_to_altloc(y) for y in  x.split("/")]) for x in g_initial})

        return g

    else:
        return nx.Graph()


def _load_connected_components(connected_components_yaml, fail_if_not_found=False):
    if fail_if_not_found and not connected_components_yaml.is_file():
        raise ValueError("File " + str(connected_components_yaml) + " does not exist")
    connected_components = {}

    if connected_components_yaml.exists():
        with open(connected_components_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic:
            for ligand_id, neighbourhood_info in dic.items():
                dtag, chain, residue, altloc, version = [dt.string_to_altloc(x) for x in ligand_id.split("+")]
                # print(f'{ligand_id} : {altloc}')
                connected_components[(dtag, chain, residue, altloc, version)] = [
                    tuple([dt.string_to_altloc(x) for x in _ligand_id.split("+")]) 
                    for _ligand_id 
                    in neighbourhood_info
                ]

    return connected_components


def _load_ligand_neighbourhood_transforms(ligand_neighbourhood_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not ligand_neighbourhood_transforms_yaml.is_file():
        raise ValueError("File " + str(ligand_neighbourhood_transforms_yaml) + " does not exist")
    ligand_neighbourhood_transforms = {}
    if ligand_neighbourhood_transforms_yaml.exists():
        with open(ligand_neighbourhood_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for ligand_transform_key, ligand_transform in dic.items():
            ligand_1_id, ligand_2_id = ligand_transform_key.split("~")
            dtag_1, chain_1, residue_1, altloc_1, version_1 = [dt.string_to_altloc(x) for x in ligand_1_id.split("/")]
            dtag_2, chain_2, residue_2, altloc_2, version_2 = [dt.string_to_altloc(x) for x in ligand_2_id.split("/")]
            ligand_neighbourhood_transforms[
                ((dtag_1, chain_1, residue_1, altloc_1, version_1), (dtag_2, chain_2, residue_2, altloc_2, version_2))
            ] = dt.Transform.from_dict(ligand_transform)

    return ligand_neighbourhood_transforms


def _load_conformer_sites(conformer_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not conformer_sites_yaml.is_file():
        raise ValueError("File " + str(conformer_sites_yaml) + " does not exist")
    conformer_sites = {}
    if conformer_sites_yaml.exists():
        with open(conformer_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)
        for conformer_site_id, conformer_site_info in dic.items():
            conformer_sites[dt.string_to_altloc(conformer_site_id)] = dt.ConformerSite.from_dict(conformer_site_info)

    return conformer_sites


def _load_conformer_site_transforms(conformer_site_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not conformer_site_transforms_yaml.is_file():
        raise ValueError("File " + str(conformer_site_transforms_yaml) + " does not exist")
    conformer_site_transforms = {}
    if conformer_site_transforms_yaml.exists():
        with open(conformer_site_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for conformer_site_transform_id, conformer_site_transform_info in dic.items():
            conformer_site_1, conformer_site_2 = conformer_site_transform_id.split("~")

            conformer_site_transforms[(conformer_site_1, conformer_site_2)] = dt.Transform.from_dict(
                conformer_site_transform_info
            )

    return conformer_site_transforms


def _load_canonical_sites(canonical_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not canonical_sites_yaml.is_file():
        raise ValueError("File " + str(canonical_sites_yaml) + " does not exist")
    canonical_sites = {}
    if canonical_sites_yaml.exists():
        with open(canonical_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)

        if dic is not None:
            for canonical_site_id, canonical_site_info in dic.items():
                canonical_sites[canonical_site_id] = dt.CanonicalSite.from_dict(canonical_site_info)

    return canonical_sites


def _load_canonical_site_transforms(canonical_site_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not canonical_site_transforms_yaml.is_file():
        raise ValueError("File " + str(canonical_site_transforms_yaml) + " does not exist")
    canonical_site_transforms = {}
    if canonical_site_transforms_yaml.exists():
        with open(canonical_site_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for canonical_site_id, canonical_site_transform_info in dic.items():
            canonical_site_transforms[canonical_site_id] = dt.Transform.from_dict(canonical_site_transform_info)

    return canonical_site_transforms


def _load_xtalform_sites(xtalform_sites_yaml, fail_if_not_found=False):
    if fail_if_not_found and not xtalform_sites_yaml.is_file():
        raise ValueError("File " + str(xtalform_sites_yaml) + " does not exist")
    xtalform_sites = {}
    if xtalform_sites_yaml.exists():
        with open(xtalform_sites_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for xtalform_site_id, xtalform_site_info in dic.items():
            xtalform_sites[xtalform_site_id] = dt.XtalFormSite.from_dict(xtalform_site_info)

    return xtalform_sites


def _load_reference_stucture_transforms(reference_structure_transforms_yaml, fail_if_not_found=False):
    if fail_if_not_found and not reference_structure_transforms_yaml.is_file():
        raise ValueError("File " + str(reference_structure_transforms_yaml) + " does not exist")
    reference_structure_transforms = {}
    if reference_structure_transforms_yaml.exists():
        with open(reference_structure_transforms_yaml, "r") as f:
            dic = yaml.safe_load(f)

        for reference_structure_transform_id, reference_structure_transform_info in dic.items():
            dtag, canonical_site_id = reference_structure_transform_id.split("~")
            reference_structure_transforms[(dtag, canonical_site_id)] = dt.Transform.from_dict(
                reference_structure_transform_info
            )

    return reference_structure_transforms