def _update_graph(
    alignability_graph,
    ligand_neighbourhoods,
    ligand_neighbourhood_transforms,
):
    nodes = alignability_graph.nodes
    edges = alignability_graph.edges

    for ligand_id in ligand_neighbourhoods:
        if ligand_id not in nodes:
            alignability_graph.add_node(ligand_id)

    for to_ligand_id, from_ligand_id in ligand_neighbourhood_transforms:
        if (to_ligand_id, from_ligand_id) not in edges:
            alignability_graph.add_edge(to_ligand_id, from_ligand_id)




