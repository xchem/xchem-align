from pathlib import Path

import networkx as nx
import numpy as np
from loguru import logger

from ligand_neighbourhood_alignment.data import (
    AtomID,
    CanonicalSite,
    CanonicalSites,
    ConformerSite,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    ResidueID,
    SiteTransforms,
    Transform,
    read_canonical_sites,
    read_graph,
    read_neighbourhoods,
    read_system_data,
    save_site_transforms,
)
from ligand_neighbourhood_alignment.save_sites import save_sites
from ligand_neighbourhood_alignment.structures import get_structures, get_transform_from_residues


def get_components(g):
    cliques = list(nx.connected_components(g))
    logger.debug(f"Cliques are: {cliques}")
    return cliques


def get_residues_from_neighbourhood(n: LigandNeighbourhood):
    rids = []
    for atom in n.atoms:
        aid: AtomID = atom.atom_id
        rid = ResidueID(chain=aid.chain, residue=aid.residue)
        rids.append(rid)

    return list(set(rids))


def get_subsites_from_components(components, neighbourhoods: LigandNeighbourhoods):
    ss = []
    j = 0
    for component in components:
        rs = []
        for lid in component:
            n: LigandNeighbourhood = neighbourhoods.get_neighbourhood(lid)
            lrs: list[ResidueID] = get_residues_from_neighbourhood(n)
            rs += lrs
        s = ConformerSite(
            id=j,
            name="",
            residues=list(set(rs)),
            members=component,
            reference_ligand_id=component[0],
        )
        j += 1
        ss.append(s)

    return ss


def get_sites_from_subsites(subsites: list[ConformerSite], neighbourhoods: LigandNeighbourhoods):
    g = nx.Graph()

    # Add the nodes
    for ss in subsites:
        g.add_node(ss.id)

    # Form the site overlap matrix
    arr = np.zeros((len(subsites), len(subsites)))
    for ss1 in subsites:
        for ss2 in subsites:
            if ss1.id == ss2.id:
                continue
            v = set(ss1.residues).intersection(set(ss2.residues))
            logger.debug(f"{ss1.id} {ss2.id} {len(v)}")
            if len(v) > 5:
                arr[ss1.id, ss2.id] = 1

    # Complete the graph
    for idx, conn in np.ndenumerate(arr):
        x, y = idx

        if x == y:
            continue
        if conn:
            g.add_edge(x, y)

    logger.debug(f"Number of sites sharing residues: {len(g.edges)}")

    # Get the connected components
    cc = get_components(g)

    # Form the sites
    sites = []
    j = 0
    for component in cc:
        members = list(set(sum([subsites[j].members for j in component], start=[])))
        subsites = [subsites[j] for j in component]
        s = CanonicalSite(
            id=j,
            subsite_ids=[subsites[j].id for j in component],
            subsites=subsites,
            members=members,
            residues=list(set(sum([subsites[j].residues for j in component], start=[]))),
            reference_ligand_id=subsites[0].reference_ligand_id,
            reference_subsite_id=subsites[0].id,
            reference_subsite=subsites[0],
        )
        j += 1
        sites.append(s)

    return sites


def get_subsite_transforms(sites: CanonicalSites, structures):
    transforms = {}
    for site_id, site in zip(sites.site_ids, sites.sites):
        rss = site.reference_ligand_id.dtag
        rs = site.residues
        srs = structures[rss]

        for ssid, ss in zip(site.subsite_ids, site.subsites):
            ssr = ss.reference_ligand_id.dtag
            ssrs = structures[ssr]
            transform = get_transform_from_residues(rs, srs, ssrs)
            transforms[(site_id, 0, ssid)] = Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist())

    return transforms


def get_site_transforms(sites: CanonicalSites, structures):
    transforms = {}
    rs = sites.reference_site
    rsid = sites.reference_site_id

    rss = structures[rs.reference_ligand_id.dtag]
    ref_site_all_ress = [
        ResidueID(chain=chain.name, residue=res.seqid.num) for model in rss for chain in model for res in chain
    ]

    for site_id, site in zip(sites.site_ids, sites.sites):
        srs = site.reference_ligand_id.dtag
        site_structure = structures[srs]

        transform = get_transform_from_residues(ref_site_all_ress, rss, site_structure)
        transforms[(rsid, site_id)] = Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist())

    return transforms


def update_subsites_from_components(sites, connected_components, neighbourhoods):
    ...


def update_sites_from_subsites(sites, updated_subsites, neighbourhoods):
    ...


def _update_sites(_source_dir: Path):
    logger.info(f"Source dir: {_source_dir}")
    g = read_graph(_source_dir)
    neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
    logger.info(f"Number of neighbourhoods: {len(neighbourhoods.ligand_ids)}")

    # Get the current sites
    sites: CanonicalSites = read_canonical_sites(_source_dir)

    # Get the (potentially new) data
    system_data = read_system_data(_source_dir)

    # Get the connected components
    logger.info("Getiting connected components...")
    connected_components = get_components(g)
    logger.info(f"Number of connected components: {len(connected_components)}")

    # Update the subsites with new component members
    # Create new subsites if necessary
    logger.info("Geting sites...")
    updated_subsites: list[ConformerSite] = update_subsites_from_components(
        sites, connected_components, neighbourhoods
    )
    logger.info(f"Number of subsites: {len(updated_subsites)}")

    # Update sites with new subsites
    logger.info("Getting sites...")
    updated_sites: CanonicalSites = update_sites_from_subsites(sites, updated_subsites, neighbourhoods)
    logger.info(f"Number of sites: {len(updated_sites.sites)}")

    save_sites(updated_sites, _source_dir)

    # Regenerate subsite transforms
    logger.info("Getting transfroms between subsites...")
    structures = get_structures(system_data)
    subsite_transforms = get_subsite_transforms(sites, structures)

    # Regenerate the site transforms
    logger.info("Getting transforms between sites...")
    site_transforms = get_site_transforms(sites, structures)
    site_transforms = SiteTransforms(
        canonical_site_transform_ids=[key for key in site_transforms.keys()],
        canonical_site_transforms=[tr for tr in site_transforms.values()],
        conformer_site_transform_ids=[key for key in subsite_transforms.keys()],
        conformer_site_transforms=[tr for tr in subsite_transforms.values()],
    )
    save_site_transforms(site_transforms, _source_dir)

    return sites
