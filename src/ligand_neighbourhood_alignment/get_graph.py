import networkx as nx
import numpy as np

from ligand_neighbourhood_alignment.data import LigandNeighbourhoods

# from loguru import logger
# from matplotlib import pyplot as plt


def get_graph(arr, ligand_neighbourhoods: LigandNeighbourhoods):
    ligand_ids = [ligand_id for ligand_id in ligand_neighbourhoods.ligand_ids]

    g = nx.Graph()
    for j in ligand_ids:
        g.add_node(j)

    for idx, conn in np.ndenumerate(arr):
        x, y = idx
        ligand_id_1 = ligand_ids[x]
        ligand_id_2 = ligand_ids[y]

        if x == y:
            continue
        if conn:
            g.add_edge(ligand_id_1, ligand_id_2)

    return g
