import numpy as np
import os
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import csv


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


def compute_all():
    # Laptop files
    #sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/Processed/atlas/hcpmmp1.csv"
    #lookup_txt = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"

    # Desktop files
    sc_path = "/home/nikita/Nikita_MRI/me/atlas/hcpmmp1.csv"
    lookup_txt = "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt"

    # Build lookup dictionary for node labels
    lookup = lookup_dictionary(lookup_txt)

    # Load connectivity matrix and create graph.
    sc_matrix = threshold_matrix_by_clipping(sc_path)
    G = create_graph(sc_matrix)

    # Compute all graph metrics.
    metrics = compute_metrics(G)

    # Print computed global metrics.
    print("Average clustering coefficient:", metrics['avg_clustering'])
    print("Global efficiency:", metrics['global_efficiency'])
    print("Average shortest path length:", metrics['avg_path_length'])
    print("Degree assortativity coefficient:", metrics['assortativity'])
    print("Modularity (Girvan-Newman partition):", metrics['modularity'])
    print("Number of communities detected:", len(metrics['communities']), "\n")

    # Print top nodes for several metrics.
    print_top_nodes(metrics, lookup, top_n=5)


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


def test_code():

    root = input("Paste complete file path to study folder: ")
    thr_choice = int(input("Would you to use 1) percentile clipping or 2) weight thresholding?: "))
    print("Running test code...")
    # 1) Load your lookup (optional, for labeling nodes)
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"
    #lookup_path = "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt"

    lookup = lookup_dictionary(lookup_path)

    # 2) Load your connectivity matrix
    sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"
    #sc_path = "/home/nikita/Nikita_MRI/me/atlas/hcpmmp1.csv"

    if thr_choice == 1:
        matrix = threshold_matrix_by_clipping(sc_path, 5, 99)
    elif thr_choice == 2:
        matrix, name_prefix = threshold_matrix_by_weight(sc_path, 0.2)
    # 3) Create the graph
    G = create_graph(matrix)

    # 4) Compute centrality metrics

    # 4a) Degree Centrality (non-weighted) and Strength (weighted)
    degree_centrality = nx.degree_centrality(G)
    strength = {node: sum(data['weight'] for _, _, data in G.edges(node, data=True))
                for node in G.nodes}

    # 4b) Eigenvector Centrality
    #     If your graph is large or possibly disconnected, you might need to
    #     increase max_iter or handle components separately.
    eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)

    # 4c) Betweenness Centrality
    betweenness_centrality = nx.betweenness_centrality(G, weight='weight')

    # 5) Global Reaching Centrality (based on degree by default)
    grc = global_reaching_centrality(G, centrality_func=nx.degree_centrality)

    # 6) Efficiency metrics
    global_efficiency = nx.global_efficiency(G)
    local_efficiency = nx.local_efficiency(G)

    output_dir = os.path.join(root, "group_analysis", "connectome_analysis_outputs")
    os.makedirs(output_dir, exist_ok=True)

    # 7) Save node-level metrics to a CSV file
    node_csv_file = os.path.join(output_dir, f"node_metrics_{name_prefix}.csv")
    with open(node_csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write header row
        writer.writerow(["Label", "Degree Centrality", "Strength", "Eigenvector Centrality", "Betweenness Centrality"])
        # Write metrics for each node
        for node in G.nodes():
            label = lookup.get(node, f"Node_{node}")
            d = degree_centrality.get(node, 0)
            s = strength.get(node, 0)
            e = eigenvector_centrality.get(node, 0)
            b = betweenness_centrality.get(node, 0)
            writer.writerow([label, f"{d:.4f}", f"{s:.4f}", f"{e:.4f}", f"{b:.4f}"])

    # 8) Save global metrics to a separate CSV file
    global_csv_file = os.path.join(output_dir, f"global_metrics_{name_prefix}.csv")
    with open(global_csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write header row
        writer.writerow(["Global Reaching Centrality", "Global Efficiency", "Local Efficiency"])
        # Write the metrics
        writer.writerow([f"{grc:.4f}", f"{global_efficiency:.4f}", f"{local_efficiency:.4f}"])

    print("Saved 2 CSVs")


def main():
    function_list = [
        ("Compute All", compute_all),
        ("Test Code", test_code),
    ]

    for idx, (desc, _) in enumerate(function_list, start=1):
        print(f"{idx}. {desc}")

    # Get user input for the function to execute
    choice = input("Enter the number corresponding to the function you want to execute: ")
    try:
        choice = int(choice)
        if 1 <= choice <= len(function_list):
            _, func = function_list[choice - 1]
            func()
        else:
            print("Invalid number. Please choose a valid option.")
    except ValueError:
        print("Invalid input. Please enter a number.")


if __name__ == "__main__":
    main()

