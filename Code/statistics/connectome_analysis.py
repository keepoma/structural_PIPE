import numpy as np
import os
import networkx as nx
import csv

from helpers.statistical_helpers import (lookup_dictionary, threshold_matrix_by_weight,
                                         threshold_matrix_by_clipping, create_graph,
                                         load_node_metrics_as_dataframe)
from visualization import (visualize_matrix_weights, visualize_saved_metrics,
                           plot_metric_boxplot, plot_metric_violin)


def global_reaching_centrality(graph, centrality_func=nx.degree_centrality):
    """
    Computes the Global Reaching Centrality (GRC) for a given graph,
    based on a chosen centrality measure (by default, degree centrality).

    One commonly cited definition is:
        GRC = (1 / (N-1)) * Σ (Cmax - Ci)
    where:
        Ci = centrality of node i
        Cmax = maximum centrality among all nodes
        N = total number of nodes in the graph
    """

    centrality_scores = centrality_func(graph)
    c_max = max(centrality_scores.values())
    N = len(centrality_scores)

    # Sum up (C_max - C_i) for all nodes
    grc_sum = sum((c_max - c) for c in centrality_scores.values())
    return grc_sum / (N - 1)


def compute_metrics_for_weight_threshold_range(root, sc_path,
                                               lookup_path, thresholds,
                                               binarize=False, overwrite=False):
    """
    For each threshold in 'thresholds', this function:
      1) Thresholds the connectivity matrix by weight.
      2) Creates a graph.
      3) Computes node-level and global metrics.
      4) Saves results to CSV files in an output directory.

    Returns:
      threshold_to_node_csv (dict): {threshold: path_to_node_csv, ...}
      threshold_to_global_csv (dict): {threshold: path_to_global_csv, ...}
    """

    # output directories
    output_dir = os.path.join(root, "group_analysis", "connectome_analysis_outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Load lookup (for node labels)
    lookup = lookup_dictionary(lookup_path)

    # Dictionaries to store CSV file paths
    threshold_to_node_csv = {}
    threshold_to_global_csv = {}

    # Loop over each threshold
    for wt in thresholds:
        name_prefix = str(wt).replace('.', '') + "wt"

        node_csv_file = os.path.join(output_dir, f"node_metrics_{name_prefix}.csv")
        global_csv_file = os.path.join(output_dir, f"global_metrics_{name_prefix}.csv")

        node_exists = os.path.isfile(node_csv_file)
        global_exists = os.path.isfile(global_csv_file)

        if node_exists and global_exists and not overwrite:
            print(f"[SKIP] Threshold {wt}: CSVs already exist, skipping computation.")
            threshold_to_node_csv[wt] = node_csv_file
            threshold_to_global_csv[wt] = global_csv_file
            continue

        print(f"[COMPUTE] Threshold {wt}: generating metrics...")

        matrix, _ = threshold_matrix_by_weight(sc_path, weight_threshold=wt, binarize=binarize)

        G = create_graph(matrix)

        # Compute metrics
        degree_centrality = nx.degree_centrality(G)
        strength = {
            node: sum(data['weight'] for _, _, data in G.edges(node, data=True))
            for node in G.nodes
        }
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight')
        grc = global_reaching_centrality(G, centrality_func=nx.degree_centrality)
        global_efficiency = nx.global_efficiency(G)
        local_efficiency = nx.local_efficiency(G)

        # Save node-level metrics to CSV
        with open(node_csv_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(
                ["Label", "Degree Centrality", "Strength", "Eigenvector Centrality", "Betweenness Centrality"])
            for node in G.nodes():
                label = lookup.get(node, f"Node_{node}")
                d = degree_centrality.get(node, 0)
                s = strength.get(node, 0)
                e = eigenvector_centrality.get(node, 0)
                b = betweenness_centrality.get(node, 0)
                writer.writerow([label, f"{d:.4f}", f"{s:.4f}", f"{e:.4f}", f"{b:.4f}"])

        # Save global metrics to CSV
        with open(global_csv_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Global Reaching Centrality", "Global Efficiency", "Local Efficiency"])
            writer.writerow([f"{grc:.4f}", f"{global_efficiency:.4f}", f"{local_efficiency:.4f}"])

        threshold_to_node_csv[wt] = node_csv_file
        threshold_to_global_csv[wt] = global_csv_file

        print(f"Threshold {wt}: Saved node-level metrics to {node_csv_file}")
        print(f"Threshold {wt}: Saved global metrics to {global_csv_file}\n")

    return threshold_to_node_csv, threshold_to_global_csv


def find_top_edges(sc_path, lookup_path, threshold=None, top_n=None):
    """
    Finds edges with the highest weights in the matrix.
    """

    # Load matrix and lookup
    sc = np.genfromtxt(sc_path, delimiter=',')
    lookup = lookup_dictionary(lookup_path)

    edges = []
    # Return dimension of NumPy array, so n is basically n of nodes (rows)
    n = sc.shape[0]

    """
    In a 2D NxN matrix there is a:
    - diagonal (cells i,i)
    - lower triangle (cells j < i)
    - upper triangle (cells j > i)
    """
    # Iterate over upper triangle (assuming symmetric matrix)
    # loops over every row
    for i in range(n):
        # loops only over columns greater than i
        for j in range(i + 1, n):
            w = sc[i, j]
            edges.append((i, j, w))

    # Filter by threshold (optional)
    if threshold is not None:
        edges = [edge for edge in edges if edge[2] > threshold]

    # Sort descending by weight
    edges.sort(key=lambda x: x[2], reverse=True)

    # Take top_n (optional)
    if top_n is not None and top_n < len(edges):
        edges = edges[:top_n]

    # Convert node indices to labels
    labeled_edges = []
    for (i, j, w) in edges:
        label_i = lookup.get(i, f"Node_{i}")
        label_j = lookup.get(j, f"Node_{j}")
        labeled_edges.append((label_i, label_j, w))

    return labeled_edges


def find_top_nodes_by_strength(sc_path, lookup_path, top_n=10):
    """
    Finds the nodes with the highest total connection weight.
    """

    sc = np.genfromtxt(sc_path, delimiter=',')
    lookup = lookup_dictionary(lookup_path)

    # If the matrix is NxN, the strength of node i is the sum of row i (or column i).
    strengths = np.nansum(sc, axis=1)

    # Create list of (node_index, strength)
    indexed_strengths = list(enumerate(strengths))

    # Sort descending by strength
    indexed_strengths.sort(key=lambda x: x[1], reverse=True)

    # Convert to labels
    top_nodes = []
    for i, val in indexed_strengths[:top_n]:
        label = lookup.get(i, f"Node_{i}")
        top_nodes.append((label, val))

    return top_nodes


def main():
    thresholds = np.linspace(0.1, 1, 9)

    root = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain"
    sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"

    #visualize_matrix_weights(sc_path)

    threshold_to_node_csv, threshold_to_global_csv = compute_metrics_for_weight_threshold_range(
        root=root,
        sc_path=sc_path,
        lookup_path=lookup_path,
        thresholds=thresholds,
        binarize=False,
        overwrite=False
    )
    #visualize_saved_metrics(threshold_to_node_csv, threshold_to_global_csv)

    top_edges = find_top_edges(
        sc_path=sc_path,
        lookup_path=lookup_path,
        threshold=10000,
        top_n=10
    )
    print("==== Strongest connections between: ====")
    for edge in top_edges:
        print(edge)

    top_nodes = find_top_nodes_by_strength(
        sc_path=sc_path,
        lookup_path=lookup_path,
        top_n=10
    )
    print("==== Strongest nodes: ====")
    for node in top_nodes:
        print(node)




if __name__ == "__main__":
    main()

