import numpy as np
import os
import networkx as nx
import csv
from helpers.statistical_helpers import (lookup_dictionary, threshold_matrix_by_weight,
                                         threshold_matrix_by_clipping, create_graph)
from visualization import visualize_saved_metrics


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
    # 1) Prepare output directories
    output_dir = os.path.join(root, "group_analysis", "connectome_analysis_outputs")
    os.makedirs(output_dir, exist_ok=True)

    # 2) Load lookup (for node labels)
    lookup = lookup_dictionary(lookup_path)

    # 3) Dictionaries to store CSV file paths
    threshold_to_node_csv = {}
    threshold_to_global_csv = {}

    # 4) Loop over each threshold
    for wt in thresholds:
        # Build file name prefix for CSV files
        name_prefix = str(wt).replace('.', '') + "wt"

        node_csv_file = os.path.join(output_dir, f"node_metrics_{name_prefix}.csv")
        global_csv_file = os.path.join(output_dir, f"global_metrics_{name_prefix}.csv")

        # Check if files exist
        node_exists = os.path.isfile(node_csv_file)
        global_exists = os.path.isfile(global_csv_file)

        # Decide whether to skip or compute
        if node_exists and global_exists and not overwrite:
            print(f"[SKIP] Threshold {wt}: CSVs already exist, skipping computation.")
            threshold_to_node_csv[wt] = node_csv_file
            threshold_to_global_csv[wt] = global_csv_file
            continue

        print(f"[COMPUTE] Threshold {wt}: generating metrics...")

        # 1) Threshold matrix
        matrix, _ = threshold_matrix_by_weight(sc_path, weight_threshold=wt, binarize=binarize)

        # 2) Create the graph
        G = create_graph(matrix)

        # 3) Compute metrics
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

        # 4) Save node-level metrics to CSV
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

        # 5) Save global metrics to CSV
        with open(global_csv_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Global Reaching Centrality", "Global Efficiency", "Local Efficiency"])
            writer.writerow([f"{grc:.4f}", f"{global_efficiency:.4f}", f"{local_efficiency:.4f}"])

        threshold_to_node_csv[wt] = node_csv_file
        threshold_to_global_csv[wt] = global_csv_file

        print(f"Threshold {wt}: Saved node-level metrics to {node_csv_file}")
        print(f"Threshold {wt}: Saved global metrics to {global_csv_file}\n")

    return threshold_to_node_csv, threshold_to_global_csv


def main():
    thresholds = np.linspace(0.1, 1, 9)

    root = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain"
    sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"

    threshold_to_node_csv, threshold_to_global_csv = compute_metrics_for_weight_threshold_range(
        root=root,
        sc_path=sc_path,
        lookup_path=lookup_path,
        thresholds=thresholds,
        binarize=False,
        overwrite=False
    )

    visualize_saved_metrics(threshold_to_node_csv, threshold_to_global_csv)


if __name__ == "__main__":
    main()

