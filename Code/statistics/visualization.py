import numpy as np
import os
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import pandas as pd


def visualize_matrix_weights(sc_path, bins=50):
    """
    Loads a matrix and plots a histogram of its values.
    """

    matrix = np.genfromtxt(sc_path, delimiter=',')

    matrix_c99 = matrix.copy()
    upper_bound = np.percentile(matrix_c99, 99)
    matrix_c99 = np.clip(matrix_c99, None, upper_bound)

    matrix_c95 = matrix.copy()
    upper_bound = np.percentile(matrix_c95, 95)
    matrix_c95 = np.clip(matrix_c95, None, upper_bound)

    matrix_c90 = matrix.copy()
    upper_bound = np.percentile(matrix_c90, 90)
    matrix_c90 = np.clip(matrix_c90, None, upper_bound)

    # Create side-by-side subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # Top left: non clipped
    axs[0, 0].hist(matrix, bins=bins, edgecolor='black')
    axs[0, 0].set_title("No Clipping")
    axs[0, 0].set_xlabel("Weight")
    axs[0, 0].set_ylabel("Frequency")

    # Top-right: 99% clipped
    axs[0, 1].hist(matrix_c99, bins=bins, edgecolor='black')
    axs[0, 1].set_title("99% Clipped")
    axs[0, 1].set_xlabel("Weight")
    axs[0, 1].set_ylabel("Frequency")

    # Bottom-left: 95% clipped
    axs[1, 0].hist(matrix_c95, bins=bins, edgecolor='black')
    axs[1, 0].set_title("95% Clipped")
    axs[1, 0].set_xlabel("Weight")
    axs[1, 0].set_ylabel("Frequency")

    # Bottom-right: 90% clipped
    axs[1, 1].hist(matrix_c90, bins=bins, edgecolor='black')
    axs[1, 1].set_title("90% Clipped")
    axs[1, 1].set_xlabel("Weight")
    axs[1, 1].set_ylabel("Frequency")

    plt.tight_layout()
    plt.show()


def visualize_graph(G, lookup):
    """
    Visualize the graph using a spring layout and label nodes using the lookup dictionary
    This function needs to be worked on
    """

    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(10, 8))

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=180)
    # Draw edges scaled by weight
    edges = G.edges(data=True)
    edge_widths = [data.get('weight', 1) * 1.1 for _, _, data in edges]
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray')

    # Create labels from the lookup dictionary
    labels = {node: lookup.get(node, str(node)) for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    plt.axis('off')
    plt.title("Structural Connectivity Graph")
    plt.show()


def visualize_matrix(matrix_path, clip):
    matrix = np.loadtxt(matrix_path, delimiter=",")
    title = os.path.splitext(os.path.basename(matrix_path))[0]
    if clip:
        upper_bound = np.percentile(matrix, 99)
        matrix = np.clip(matrix, None, upper_bound)
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, cmap='viridis', square=True)
    plt.title(title)
    plt.xlabel('Region Index')
    plt.ylabel('Region Index')
    plt.show()


def visualize_saved_metrics(threshold_to_node_csv, threshold_to_global_csv):
    """
    Reads the saved CSV files from 'compute_metrics_for_weight_threshold_range'
    and plots how the metrics vary with the threshold.
    - For node-level metrics, compute the average across nodes.
    - For global metrics, read the single row from each CSV file.
    """

    # Sort thresholds for consistent plotting
    sorted_thresholds = sorted(threshold_to_node_csv.keys())

    # Lists to hold aggregated data
    avg_degree_list = []
    avg_strength_list = []
    avg_eigen_list = []
    avg_betw_list = []
    grc_list = []
    global_eff_list = []
    local_eff_list = []

    # Loop over thresholds in sorted order
    for wt in sorted_thresholds:
        node_csv = threshold_to_node_csv[wt]
        global_csv = threshold_to_global_csv[wt]

        # -- Node-level metrics --
        node_df = pd.read_csv(node_csv)
        avg_degree_list.append(node_df["Degree Centrality"].mean())
        avg_strength_list.append(node_df["Strength"].mean())
        avg_eigen_list.append(node_df["Eigenvector Centrality"].mean())
        avg_betw_list.append(node_df["Betweenness Centrality"].mean())

        # -- Global metrics --
        global_df = pd.read_csv(global_csv)
        # There's only one row in each global CSV
        grc_list.append(global_df.iloc[0]["Global Reaching Centrality"])
        global_eff_list.append(global_df.iloc[0]["Global Efficiency"])
        local_eff_list.append(global_df.iloc[0]["Local Efficiency"])

    # Plot node-level metrics (excluding strength) in one figure
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_thresholds, avg_degree_list, marker='o', label="Avg Degree Centrality")
    plt.plot(sorted_thresholds, avg_eigen_list,  marker='o', label="Avg Eigenvector Centrality")
    plt.plot(sorted_thresholds, avg_betw_list,   marker='o', label="Avg Betweenness Centrality")
    plt.xlabel("Weight Threshold")
    plt.ylabel("Average Node-level Metric")
    plt.title("Node-level Metrics vs Weight Threshold (Excl. Strength)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot strength alone in another figure
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_thresholds, avg_strength_list, marker='o', color='red', label="Avg Strength")
    plt.xlabel("Weight Threshold")
    plt.ylabel("Average Strength")
    plt.title("Node Strength vs Weight Threshold")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot global metrics in a third figure
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_thresholds, grc_list, marker='o', label="Global Reaching Centrality")
    plt.plot(sorted_thresholds, global_eff_list, marker='o', label="Global Efficiency")
    plt.plot(sorted_thresholds, local_eff_list, marker='o', label="Local Efficiency")
    plt.xlabel("Weight Threshold")
    plt.ylabel("Global Metric Value")
    plt.title("Global Metrics vs Weight Threshold")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_metric_boxplot(df, metric="Degree Centrality"):
    """
    Plots a boxplot of the given 'metric' against the 'Threshold' column.
    Expects that 'df' has columns: ['Threshold', metric].
    """
    # Convert threshold to string so Seaborn treats it as a categorical variable
    df["Threshold"] = df["Threshold"].astype(str)

    plt.figure(figsize=(8, 6))
    sns.boxplot(x="Threshold", y=metric, data=df)
    plt.title(f"Distribution of {metric} by Threshold (Boxplot)")
    plt.xlabel("Threshold")
    plt.ylabel(metric)
    plt.tight_layout()
    plt.show()


def plot_metric_violin(df, metric="Degree Centrality"):
    """
    Plots a violin plot of the given 'metric' against the 'Threshold' column.
    """
    # Convert threshold to string so Seaborn treats it as a categorical variable
    df["Threshold"] = df["Threshold"].astype(str)

    plt.figure(figsize=(8, 6))
    sns.violinplot(x="Threshold", y=metric, data=df, inner="box")
    plt.title(f"Distribution of {metric} by Threshold (Violin)")
    plt.xlabel("Threshold")
    plt.ylabel(metric)
    plt.tight_layout()
    plt.show()