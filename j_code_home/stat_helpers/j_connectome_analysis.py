import numpy as np
import os
import networkx as nx
import csv
from scipy import stats
from stat_helpers.j_statistical_helpers import (lookup_dictionary, threshold_matrix_by_weight,
                                         threshold_matrix_by_clipping, create_graph,
                                         load_node_metrics_as_dataframe)
from connectome_visualization import (visualize_matrix_weights, visualize_saved_metrics,
                                      plot_metric_boxplot, plot_metric_violin,
                                      visualize_matrix_comparison, visualize_matrix,
                                      visualize_matrix_side_by_side)





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



def compute_connectivity_metrics(matrix, csv_path=None):
    """
    Reads a structural connectivity matrix from a file and computes several metrics.
    Returns:
        dict: A dictionary with keys 'mean', 'median', 'q1', 'q3', 'std', and 'variance'.
              (You can extend it to include mode or other statistics as needed.)

    If csv_path is provided, the metrics are also saved to the specified CSV file.
    """

    # Flatten the matrix into a 1D array for computing overall statistics
    flat_matrix = matrix.flatten()

    mean_val = np.mean(flat_matrix)
    median_val = np.median(flat_matrix)

    # Compute quartiles: 25th (Q1) and 75th (Q3) percentiles
    q1 = np.percentile(flat_matrix, 25)
    q3 = np.percentile(flat_matrix, 75)

    # Standard deviation and variance
    std_val = np.std(flat_matrix)
    var_val = np.var(flat_matrix)

    metrics = {
        'mean': mean_val,
        'median': median_val,
        'q1': q1,
        'q3': q3,
        'std': std_val,
        'variance': var_val
    }

    # Optionally save the metrics to a CSV file
    if csv_path is not None:
        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write a header row (optional)
            writer.writerow(["metric", "value"])
            for key, value in metrics.items():
                writer.writerow([key, value])

    return metrics


def compute_metrics_for_weight_threshold_range(paths, sc_path,
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
    con_stats_dir = paths["con_stats_dir"]
    os.makedirs(paths["con_stats_dir"], exist_ok=True)


    # Load lookup (for node labels)
    lookup = lookup_dictionary(lookup_path)

    # Dictionaries to store CSV file paths
    threshold_to_node_csv = {}
    threshold_to_global_csv = {}

    # Loop over each threshold
    for wt in thresholds:
        name_prefix = str(wt).replace('.', '') + "wt"

        node_csv_file = os.path.join(con_stats_dir, f"node_metrics_{name_prefix}.csv")
        global_csv_file = os.path.join(con_stats_dir, f"global_metrics_{name_prefix}.csv")

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


def find_top_nodes_by_strength(matrix, lookup_path, top_n=None, csv_path=None):
    """
    Finds the nodes sorted by total connection weight.
    If top_n is provided, returns only the top_n nodes; otherwise, returns all nodes.
    If csv_path is provided, saves the results to the specified CSV file.
    """
    # Calculate node strengths
    strengths = np.nansum(matrix, axis=1)

    # Load your lookup dictionary
    lookup = lookup_dictionary(lookup_path)

    # Create (index, strength) pairs
    indexed_strengths = list(enumerate(strengths))

    # Sort descending by strength
    indexed_strengths.sort(key=lambda x: x[1], reverse=True)

    # Select the top_n nodes (or all if top_n is None)
    selected = indexed_strengths if top_n is None else indexed_strengths[:top_n]

    # Convert indices to labels
    top_nodes = [(lookup.get(i, f"Node_{i}"), val) for i, val in selected]

    # If a csv_path is provided, write the results to a CSV file
    if csv_path is not None:
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Optional header
            writer.writerow(["Node", "Strength"])
            for label, val in top_nodes:
                writer.writerow([label, val])

    return top_nodes
