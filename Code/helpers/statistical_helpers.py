import numpy as np
import networkx as nx


def lookup_dictionary(lookup_txt):
    """
    Build a lookup dictionary mapping node IDs to label names
    """

    lookup = {}
    with open(lookup_txt, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                node_id = int(parts[0])
                label = parts[1]
                lookup[node_id - 1] = label  # Adjust to 0-indexed
            except Exception as e:
                print("Skipping line due to error:", line, e)
    return lookup


def threshold_matrix_by_clipping(sc_path, lower_pct, upper_pct, binarize=False):
    """
    Load the structural connectivity matrix from a CSV file and clip values
    """

    sc = np.genfromtxt(sc_path, delimiter=',')

    # 1) Clip (not remove) extreme high values at upper_pct percentile
    high_val = np.percentile(sc, upper_pct)
    sc[sc > high_val] = high_val

    # 2) Threshold at the lower_pct percentile
    low_val = np.percentile(sc, lower_pct)
    sc[sc < low_val] = 0

    # (Optional) Binarize
    if binarize:
        sc[sc > 0] = 1

    return sc

def threshold_matrix_by_weight(sc_path, weight_threshold=0.2, binarize=False):
    """
    Zeros out all edge weights below the specified weight_threshold.
    """

    sc = np.genfromtxt(sc_path, delimiter=',')
    sc[sc < weight_threshold] = 0.0
    if binarize:
        sc[sc > 0] = 1.0

    name_prefix = str(weight_threshold).replace('.', '') + "wt"
    return sc, name_prefix


def create_graph(matrix):
    """
    Create a graph from the connectivity matrix.
    Undirected and weighted from current code
    """

    G = nx.from_numpy_array(matrix)

    # Remove self-loops if present
    G.remove_edges_from(nx.selfloop_edges(G))
    return G