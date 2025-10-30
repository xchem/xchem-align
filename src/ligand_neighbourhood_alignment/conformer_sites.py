import networkx as nx

from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment import alignment_heirarchy


def _get_connected_components(alignability_graph, clusters, max_path_length=2):
    """
    Construct neighbourhoods around the most connected neighbourhoods by some max path length,


    """

    # Get the graph of short paths
    path = dict(nx.all_pairs_shortest_path(alignability_graph))
    path_lengths = {(source, target): len(path[source][target]) for source in path for target in path[source]}
    H = nx.Graph()
    for node in path:
        H.add_node(node)
    for source, target in path_lengths:
        if path_lengths[(source, target)] <= max_path_length:
            H.add_edge(source, target)

    #
    degrees = dict(nx.degree(H))

    # Replay cluster cores
    used = [member for cluster in clusters for member in cluster]
    for x in clusters:
        for target in H.nodes:
            if target in used:
                continue
            if (x, target) not in path_lengths:
                continue
            if path_lengths[(x, target)] <= 2:
                used.append(target)
                clusters[x].append(target)

    # Now go through any new ligands that are not yet connected, constructing clusters for them
    for x in sorted(degrees, key=lambda _x: degrees[_x], reverse=True):
        if x in used:
            continue
        clusters[x] = []
        # print(f"f{x} : {degrees[x]}")

        # for n in G.neighbors(x):
        # used.append(n)
        for target in H.nodes:
            if target in used:
                continue
            if (x, target) not in path_lengths:
                continue
            if path_lengths[(x, target)] <= 2:
                used.append(target)
                clusters[x].append(target)

    return clusters




def _update_conformer_sites(
    conformer_sites: dict[str, dt.ConformerSite],
    connected_component_id: tuple[str, str, str, str],
    connected_component: list[tuple[str, str, str, str]],
    neighbourhoods: dict[tuple[str, str, str, str], dt.Neighbourhood],
    structures,
    assemblies,
    xtalforms,
    xtalform_assignments,
):
    matched = False
    # Check each old conformer site for overlap in membership, and if so update its members
    for conformer_site_id, conformer_site in conformer_sites.items():
        num_overlaps = set(connected_component).intersection(set(conformer_site.members))
        if len(num_overlaps) > 0:
            matched = True
            # Match, add each new ligand id to the conformer site's members
            for lid in connected_component:
                if lid not in conformer_site.members:
                    conformer_site.members.append(lid)

    # Otherwise create a new conformer site
    if not matched:
        residues = []
        residues_aligned = []
        for lid in connected_component:
            st = structures[lid[0]]
            if lid != connected_component_id:
                continue
            for atom_id in neighbourhoods[lid].atoms:
                biochain = alignment_heirarchy._chain_to_biochain(
                    atom_id[0], xtalforms[xtalform_assignments[lid[0]]], assemblies
                )
                residues.append((atom_id[0], atom_id[1], st[0][atom_id[0]][atom_id[1]][0].name))
                residues_aligned.append((biochain, atom_id[1], st[0][atom_id[0]][atom_id[1]][0].name))
        conformer_site = dt.ConformerSite(
            [x for x in set(residues)],
            [x for x in set(residues_aligned)],
            connected_component,
            # [x for x in connected_component][0]
            connected_component_id,
        )
        conformer_site_id = "+".join(
            [dt.altloc_to_string(x) for x in conformer_site.reference_ligand_id]
            )
        conformer_sites[conformer_site_id] = conformer_site


