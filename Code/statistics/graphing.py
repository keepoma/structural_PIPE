import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from connectome_analysis import (threshold_matrix_by_weight, threshold_matrix_by_clipping,
                                 create_graph, lookup_dictionary,
                                 global_reaching_centrality)


def visualize_weight_thresholding_impact():
    """
    Visualize how varying the weight threshold impacts the network measures.
    Uses the same file paths and metric computations as in your test_code() function.
    """
    # Use the same lookup and connectivity matrix paths as in test_code()
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"
    sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"

    # Load lookup dictionary
    lookup = lookup_dictionary(lookup_path)

    # Define the range of weight thresholds (for example, from 0.1 to 0.5)
    weight_thresholds = np.linspace(0.1, 1, 18)

    # Prepare lists to collect metrics
    avg_degree_list = []
    avg_strength_list = []
    avg_eigen_list = []
    avg_betw_list = []
    grc_list = []
    global_eff_list = []
    local_eff_list = []

    for wt in weight_thresholds:
        # Threshold the matrix by weight (this returns the modified matrix and a name prefix)
        matrix, _ = threshold_matrix_by_weight(sc_path, weight_threshold=wt)
        # Create graph
        G = create_graph(matrix)

        # Compute node-level metrics:
        degree_centrality = nx.degree_centrality(G)
        strength = {node: sum(data['weight'] for _, _, data in G.edges(node, data=True))
                    for node in G.nodes}
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight')

        # Compute global metrics:
        grc_val = global_reaching_centrality(G, centrality_func=nx.degree_centrality)
        global_eff = nx.global_efficiency(G)
        local_eff = nx.local_efficiency(G)

        # Average the node-level metrics over all nodes
        avg_degree = np.mean(list(degree_centrality.values()))
        avg_str = np.mean(list(strength.values()))
        avg_eigen = np.mean(list(eigenvector_centrality.values()))
        avg_betw = np.mean(list(betweenness_centrality.values()))

        avg_degree_list.append(avg_degree)
        avg_strength_list.append(avg_str)
        avg_eigen_list.append(avg_eigen)
        avg_betw_list.append(avg_betw)
        grc_list.append(grc_val)
        global_eff_list.append(global_eff)
        local_eff_list.append(local_eff)

    # Plot the results: two panels for node-level and global metrics
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))

    # Node-level metrics plot
    axs[0].plot(weight_thresholds, avg_degree_list, marker='o', label="Avg Degree Centrality")
    axs[0].plot(weight_thresholds, avg_strength_list, marker='o', label="Avg Strength")
    axs[0].plot(weight_thresholds, avg_eigen_list, marker='o', label="Avg Eigenvector Centrality")
    axs[0].plot(weight_thresholds, avg_betw_list, marker='o', label="Avg Betweenness Centrality")
    axs[0].set_xlabel("Weight Threshold")
    axs[0].set_ylabel("Average Node-level Metric")
    axs[0].set_title("Node-level Metrics vs Weight Threshold")
    axs[0].legend()

    # Global metrics plot
    axs[1].plot(weight_thresholds, grc_list, marker='o', label="Global Reaching Centrality")
    axs[1].plot(weight_thresholds, global_eff_list, marker='o', label="Global Efficiency")
    axs[1].plot(weight_thresholds, local_eff_list, marker='o', label="Local Efficiency")
    axs[1].set_xlabel("Weight Threshold")
    axs[1].set_ylabel("Global Metric Value")
    axs[1].set_title("Global Metrics vs Weight Threshold")
    axs[1].legend()

    plt.tight_layout()
    plt.show()


def visualize_clipping_impact():
    """
    Visualize how varying the lower clipping percentile impacts the network measures.
    Uses the same file paths and metric computations as in your test_code() function.
    """
    # Use the same lookup and connectivity matrix paths as in test_code()
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"
    sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"

    # Load lookup dictionary
    lookup = lookup_dictionary(lookup_path)

    # Define the range for the lower clipping percentile (for example, 1% to 10%)
    lower_percentiles = np.linspace(1, 10, 10)
    fixed_upper_pct = 99  # Using the same fixed upper percentile as in test_code()

    # Prepare lists to collect metrics
    avg_degree_list = []
    avg_strength_list = []
    avg_eigen_list = []
    avg_betw_list = []
    grc_list = []
    global_eff_list = []
    local_eff_list = []

    for lower_pct in lower_percentiles:
        # Threshold the matrix by clipping using the current lower percentile
        matrix = threshold_matrix_by_clipping(sc_path, lower_pct=lower_pct, upper_pct=fixed_upper_pct)
        # Create graph
        G = create_graph(matrix)

        # Compute node-level metrics:
        degree_centrality = nx.degree_centrality(G)
        strength = {node: sum(data['weight'] for _, _, data in G.edges(node, data=True))
                    for node in G.nodes}
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight')

        # Compute global metrics:
        grc_val = global_reaching_centrality(G, centrality_func=nx.degree_centrality)
        global_eff = nx.global_efficiency(G)
        local_eff = nx.local_efficiency(G)

        # Average the node-level metrics over all nodes
        avg_degree = np.mean(list(degree_centrality.values()))
        avg_str = np.mean(list(strength.values()))
        avg_eigen = np.mean(list(eigenvector_centrality.values()))
        avg_betw = np.mean(list(betweenness_centrality.values()))

        avg_degree_list.append(avg_degree)
        avg_strength_list.append(avg_str)
        avg_eigen_list.append(avg_eigen)
        avg_betw_list.append(avg_betw)
        grc_list.append(grc_val)
        global_eff_list.append(global_eff)
        local_eff_list.append(local_eff)

    # Plot the results: two panels for node-level and global metrics
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))

    # Node-level metrics plot
    axs[0].plot(lower_percentiles, avg_degree_list, marker='o', label="Avg Degree Centrality")
    axs[0].plot(lower_percentiles, avg_strength_list, marker='o', label="Avg Strength")
    axs[0].plot(lower_percentiles, avg_eigen_list, marker='o', label="Avg Eigenvector Centrality")
    axs[0].plot(lower_percentiles, avg_betw_list, marker='o', label="Avg Betweenness Centrality")
    axs[0].set_xlabel("Lower Clipping Percentile")
    axs[0].set_ylabel("Average Node-level Metric")
    axs[0].set_title("Node-level Metrics vs Lower Clipping Percentile")
    axs[0].legend()

    # Global metrics plot
    axs[1].plot(lower_percentiles, grc_list, marker='o', label="Global Reaching Centrality")
    axs[1].plot(lower_percentiles, global_eff_list, marker='o', label="Global Efficiency")
    axs[1].plot(lower_percentiles, local_eff_list, marker='o', label="Local Efficiency")
    axs[1].set_xlabel("Lower Clipping Percentile")
    axs[1].set_ylabel("Global Metric Value")
    axs[1].set_title("Global Metrics vs Lower Clipping Percentile")
    axs[1].legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    visualize_weight_thresholding_impact()
