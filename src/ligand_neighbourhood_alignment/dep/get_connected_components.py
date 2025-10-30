import networkx as nx
import numpy as np
from loguru import logger
from matplotlib import pyplot as plt


def get_connected_components(arr, ligand_neighbourhoods):
    ligand_ids = [ligand_id for ligand_id in ligand_neighbourhoods]

    g = nx.Graph()
    for x in ligand_ids:
        g.add_node(x)

    for idx, conn in np.ndenumerate(arr):
        x, y = idx
        ligand_id_1 = ligand_ids[x]
        ligand_id_2 = ligand_ids[y]

        if x == y:
            continue
        if conn:
            g.add_edge(ligand_id_1, ligand_id_2)

    cliques = list(nx.find_cliques(g))
    logger.debug(f"Cliques are: {cliques}")

    fig, ax = plt.subplots()

    nx.draw(g, ax=ax)
    fig.savefig("./graph.png")

    # logger.debug(g)

    return cliques


def get_connected_components_connected(arr, ligand_neighbourhoods):
    ligand_ids = [ligand_id for ligand_id in ligand_neighbourhoods]

    g = nx.Graph()
    for x in ligand_ids:
        g.add_node(x)

    for idx, conn in np.ndenumerate(arr):
        x, y = idx
        ligand_id_1 = ligand_ids[x]
        ligand_id_2 = ligand_ids[y]

        if x == y:
            continue
        if conn:
            g.add_edge(ligand_id_1, ligand_id_2)

    cliques = list(nx.connected_components(g))
    logger.debug(f"Cliques are: {cliques}")

    fig, ax = plt.subplots()

    nx.draw(g, ax=ax)
    fig.savefig("./graph.png")

    return cliques, g
