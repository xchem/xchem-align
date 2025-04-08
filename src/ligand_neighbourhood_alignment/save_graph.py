from pathlib import Path

import networkx as nx


def save_graph(g, path: Path):
    nx.write_gml(g, str(path), stringizer=lambda x: x.to_string())
