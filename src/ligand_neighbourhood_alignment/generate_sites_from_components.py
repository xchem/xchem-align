from pathlib import Path

import networkx as nx
import numpy as np
from loguru import logger

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment import alignment_heirarchy
from ligand_neighbourhood_alignment.data import (
    AssignedXtalForms,
    AtomID,
    CanonicalSite,
    CanonicalSites,
    ConformerSite,
    ConformerSites,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    ResidueID,
    SiteTransforms,
    Transform,
    XtalForms,
    XtalFormSite,
    XtalFormSites,
    read_graph,
    read_neighbourhoods,
    read_system_data,
    save_site_transforms,
)

# from ligand_neighbourhood_alignment.save_sites import save_sites
from ligand_neighbourhood_alignment.structures import (
    get_structures,
    get_transform_from_residues,
    _get_transform_from_residues,
)


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


def get_conformer_sites_from_components(components, neighbourhoods: LigandNeighbourhoods):
    ss = []
    j = 0
    for component in components:
        rs = []
        for lid in component:
            n: LigandNeighbourhood = neighbourhoods.get_neighbourhood(lid)
            lrs: list[ResidueID] = get_residues_from_neighbourhood(n)
            rs += lrs
        component_members = list(component)
        s = ConformerSite(
            id=j,
            name="",
            residues=list(set(rs)),
            members=component_members,
            reference_ligand_id=component_members[0],
        )
        j += 1
        ss.append(s)

    return ConformerSites(conformer_sites={site.id: site for site in ss})


def get_sites_from_conformer_sites(conformer_sites: ConformerSites, neighbourhoods: LigandNeighbourhoods):
    g = nx.Graph()

    # Add the nodes
    for ss_id, ss in conformer_sites.iter():
        g.add_node(ss.id)

    conformer_sites_dict = {ss.id: ss for ss_id, ss in conformer_sites.iter()}

    # Form the site overlap matrix
    arr = np.zeros(
        (
            len(conformer_sites.conformer_sites),
            len(conformer_sites.conformer_sites),
        )
    )
    for ss_id1, ss1 in conformer_sites.iter():
        for ss_id2, ss2 in conformer_sites.iter():
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
    logger.debug((f"Conformer Sites: {conformer_sites}" f" {len(conformer_sites.conformer_sites)}"))

    for component in cc:
        logger.debug(f"Component: {component} {len(component)}")
        _members = list(
            set(
                sum(
                    [conformer_sites_dict[k].members for k in component],
                    start=[],
                )
            )
        )
        _subsites = [conformer_sites_dict[k] for k in component]

        s = CanonicalSite(
            id=j,
            subsite_ids=[conformer_sites_dict[k].id for k in component],
            subsites=_subsites,
            members=_members,
            residues=list(
                set(
                    sum(
                        [conformer_sites.conformer_sites[j].residues for j in component],
                        start=[],
                    )
                )
            ),
            reference_ligand_id=conformer_sites.conformer_sites[0].reference_ligand_id,
            reference_subsite_id=conformer_sites.conformer_sites[0].id,
            reference_subsite=conformer_sites.conformer_sites[0],
        )
        logger.debug(f"Canonical site: {j} has {len(s.subsites)} conformer sites")
        j += 1
        sites.append(s)

    return sites


def get_xtalform_sites_from_canonical_sites(
    canonical_sites: CanonicalSites,
    assigned_xtalforms: AssignedXtalForms,
    xtalforms: XtalForms,
):
    """
    Each canonical site may occur in several forms, depending on the
    combination of real and artefact chains involved. This function
    partitions the site into xtalform sites, each of which contains
    a single possible crystallographic context i.e. each
    crystallographicly non-identical chain needs its own site.

    Binding may occur to multiple crystallographically non-
    identical chains in the same crystalform, giving rise to
    multiple xtalform sites for that dataset. The easiest way to
    check whether two chains are identical is whether one is
    related to the other by a transform in the crystal form
    assembly.
    """

    xtalform_site_num = 0
    xtalform_sites: dict[tuple, XtalFormSite] = {}

    for site_id, site in canonical_sites.iter():
        # site_residues = site.residues
        for ligand_id in site.members:
            chain = ligand_id.chain
            dtag = ligand_id.dtag
            xtalform_id = assigned_xtalforms.get_xtalform_id(dtag)
            xtalform = xtalforms.get_xtalform(xtalform_id)

            # Determine which crystallographic chain the ligand is part of
            # by finding the chain that generated it (normally the same chain)
            for assembly_id, assembly in xtalform.assemblies.items():
                for generator_id, generator in assembly.generators.items():
                    if chain == generator.reference_chain:
                        crystallographic_chain = generator.chain

            # Define the site key
            xtalform_site_key = (site_id, xtalform_id, crystallographic_chain)

            # Check if the xtalform assembly pair has a site, add if so
            if xtalform_site_key in xtalform_sites:
                if ligand_id not in xtalform_sites[xtalform_site_key].members:
                    xtalform_sites[xtalform_site_key].members.append(ligand_id)

            # Create a new xtalform site for the xtalform-assembly pair
            else:
                xtalform_sites[xtalform_site_key] = XtalFormSite(
                    id=xtalform_site_num,
                    site_id=site_id,
                    xtalform_id=xtalform_id,
                    crystallographic_chain=crystallographic_chain,
                    members=[
                        ligand_id,
                    ],
                )
                xtalform_site_num += 1

    return XtalFormSites(xtalform_sites={xtalform_site.id: xtalform_site for xtalform_site in xtalform_sites.values()})


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


def _update_conformer_site_transforms(
    conformer_site_transforms,
    canonical_site: dt.CanonicalSite,
    conformer_sites: dict[str, dt.ConformerSite],
    structures,
):
    ref_conformer_site = conformer_sites[canonical_site.reference_conformer_site_id]
    ref_conformer_site_residues = ref_conformer_site.residues

    for conformer_site_id in canonical_site.conformer_site_ids:
        key = (canonical_site.reference_conformer_site_id, conformer_site_id)
        if key not in conformer_site_transforms:
            conformer_site = conformer_sites[conformer_site_id]
            # conformer_site_residues = conformer_site.residues

            transform = _get_transform_from_residues(
                [(x[0], x[1]) for x in canonical_site.residues],
                structures[conformer_site.reference_ligand_id[0]],
                structures[ref_conformer_site.reference_ligand_id[0]],
            )

            conformer_site_transforms[key] = dt.Transform(transform.vec.tolist(), transform.mat.tolist())

    # transforms = {}
    # for site_id, site in zip(sites.site_ids, sites.sites):
    #     rss = site.reference_ligand_id.dtag
    #     rs = site.residues
    #     srs = structures[rss]
    #
    #     for ssid, ss in zip(site.subsite_ids, site.subsites):
    #         ssr = ss.reference_ligand_id.dtag
    #         ssrs = structures[ssr]
    #         transform = get_transform_from_residues(rs, srs, ssrs)
    #         transforms[(site_id, 0, ssid)] = Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist())

    # return transforms


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


def _update_canonical_site_transforms(
    canonical_site_transforms: dict[str, dt.Transform],
    canonical_site_id,
    canonical_site: dt.CanonicalSite,
    # canonical_sites: dict[str, dt.CanonicalSite],
    conformer_sites: dict[str, dt.ConformerSite],
    structures,
):
    rss = structures[canonical_site.global_reference_dtag]
    ref_site_all_ress = [(chain.name, res.seqid.num) for model in rss for chain in model for res in chain]

    srs = conformer_sites[canonical_site.reference_conformer_site_id].reference_ligand_id[0]
    site_structure = structures[srs]

    transform = _get_transform_from_residues(ref_site_all_ress, rss, site_structure)
    canonical_site_transforms[canonical_site_id] = dt.Transform(
        transform.vec.tolist(),
        transform.mat.tolist(),
    )


def _update_reference_structure_transforms(
    reference_structure_transforms,
    key,
    structures,
    canonical_site: dt.CanonicalSite,
    conformer_sites: dict[str, dt.ConformerSite],
    assemblies,
    xtalforms,
    dataset_assignments,
    xtalform_sites,
    canonical_site_id,
):
    # # Get the biochain of the canonical site
    site_reference_ligand_id = conformer_sites[canonical_site.reference_conformer_site_id].reference_ligand_id
    site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
    site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]
    for xsid, _xtalform_site in xtalform_sites.items():
        _xtalform_id = _xtalform_site.xtalform_id
        if (
            (_xtalform_id == site_reference_ligand_xtalform_id)
            & (_xtalform_site.canonical_site_id == canonical_site_id)
            & (site_reference_ligand_id in _xtalform_site.members)
        ):
            xtalform_site = _xtalform_site
    site_chain = xtalform_site.crystallographic_chain
    canonical_site_biochain = alignment_heirarchy._chain_to_biochain(
        site_chain, site_reference_ligand_xtalform, assemblies
    )

    # Determine whether the biochain is shared, and if not skip
    reference_structure = structures[key[0]]
    reference_structure_xtalform = xtalforms[dataset_assignments[key[0]]]
    xtalform_chains = [
        chain for assembly in reference_structure_xtalform.assemblies.values() for chain in assembly.chains
    ]
    reference_structure_biochains = {
        chain: alignment_heirarchy._chain_to_biochain(chain, reference_structure_xtalform, assemblies)
        for chain in xtalform_chains
    }
    reference_structure_biochains_inv = {v: k for k, v in reference_structure_biochains.items()}

    # # Get the reference residues whose image is the canonical site biochain aligned residues
    # alignment_residues_ref_st = []
    # alignment_residues_mov_st = []
    # for rid in canonical_site.residues:
    #     chain, res = rid[0], rid[1]
    #     biochain = alignment_heirarchy._chain_to_biochain(
    #         chain,
    #         site_reference_ligand_xtalform,
    #         assemblies
    #     )
    #     reference_structure_chain = reference_structure_biochains_inv[biochain]
    #     alignment_residues_ref_st.append((reference_structure_chain, res))
    #     alignment_residues_mov_st.append((reference_structure_chain, res))

    # # Align the reference to the biochain reference using the canonical site residues
    alignment_residues_ref_st = []
    alignment_residues_mov_st = []
    core_chain = reference_structure_biochains_inv[canonical_site_biochain]
    canonical_site_residues = canonical_site.residues
    for rid in canonical_site_residues:
        chain, res = rid[0], rid[1]
        biochain = alignment_heirarchy._chain_to_biochain(chain, site_reference_ligand_xtalform, assemblies)
        if biochain != canonical_site_biochain:
            continue
        alignment_residues_ref_st.append((chain, res))
        alignment_residues_mov_st.append((core_chain, res))
    transform = _get_transform_from_residues(
        alignment_residues_ref_st,
        structures[site_reference_ligand_id[0]],
        reference_structure,
        other_rs=alignment_residues_mov_st,
    )

    reference_structure_transforms[key] = dt.Transform(
        transform.vec.tolist(),
        transform.mat.tolist(),
    )
    return transform

    # DEP
    # ress = [(x[0], x[1]) for x in canonical_site.residues]
    # to_structure = structures[conformer_sites[canonical_site.reference_conformer_site_id].reference_ligand_id[0]]
    # from_structure = structures[key[0]]
    # transform = _get_transform_from_residues(ress, to_structure, from_structure, other_rs=other_rs)
    # reference_structure_transforms[key] = dt.Transform(
    #     transform.vec.tolist(),
    #     transform.mat.tolist(),
    # )
    # return transform


def _generate_sites_from_components(_source_dir: Path):
    logger.info(f"Source dir: {_source_dir}")
    g = read_graph(_source_dir)
    neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
    logger.info(f"Number of neighbourhoods: {len(neighbourhoods.ligand_ids)}")

    system_data = read_system_data(_source_dir)

    # Get the assemblies
    # assemblies = Assemblies.read(_source_dir / constants.ASSEMBLIES_FILE_NAME)
    assigned_xtalforms = AssignedXtalForms.read(_source_dir / constants.ASSIGNED_XTALFORMS_FILE_NAME)

    # Get the xtalforms
    xtalforms = XtalForms.read(_source_dir / constants.XTALFORMS_FILE_NAME)

    # Get the connected components
    logger.info("Getiting connected components...")
    connected_components = get_components(g)
    logger.info(f"Number of connected components: {len(connected_components)}")

    # Get the subsites from the connected components with overlap
    logger.info("Geting sites...")
    conformer_sites: ConformerSites = get_conformer_sites_from_components(connected_components, neighbourhoods)
    logger.info(f"Number of subsites: {len(conformer_sites.conformer_sites)}")

    # Save the conformer sites
    conformer_sites.save(_source_dir / constants.CONFORMER_SITE_FILE)

    # Merge the connected components with shared residues into sites
    logger.info("Getting sites...")
    _sites: list[CanonicalSite] = get_sites_from_conformer_sites(conformer_sites, neighbourhoods)
    logger.info(f"Number of sites: {len(_sites)}")

    canonical_sites: CanonicalSites = CanonicalSites(
        site_ids=[s.id for s in _sites],
        sites=_sites,
        reference_site=_sites[0],
        reference_site_id=_sites[0].id,
    )
    canonical_sites.save(_source_dir / constants.CANONICAL_SITE_FILE)

    # Get the xtalform sites
    xtalform_sites: XtalFormSites = get_xtalform_sites_from_canonical_sites(
        canonical_sites,
        assigned_xtalforms,
        xtalforms,
        # assemblies,
    )
    xtalform_sites.save(_source_dir / constants.XTALFORM_SITE_FILE)

    # Get the subsite transforms
    logger.info("Getting transfroms between subsites...")
    structures = get_structures(system_data)
    subsite_transforms = get_subsite_transforms(canonical_sites, structures)

    # Get the site transforms
    logger.info("Getting transforms between sites...")
    site_transforms = get_site_transforms(canonical_sites, structures)
    site_transforms = SiteTransforms(
        canonical_site_transform_ids=list(site_transforms.keys()),
        canonical_site_transforms=list(site_transforms.values()),
        conformer_site_transform_ids=list(subsite_transforms.keys()),
        conformer_site_transforms=list(subsite_transforms.values()),
    )
    save_site_transforms(site_transforms, _source_dir)

    return canonical_sites
