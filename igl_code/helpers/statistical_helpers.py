import os
from Code.helpers.helpers import get_subject_dirs
import numpy as np
import networkx as nx
import pandas as pd


def chose_workspace(workspace):
    # Getting custom user input
    if workspace == 1:
        root = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain"
        sc_path = "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/atlas/hcpmmp1.csv"
        lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/hcpmmp1_ordered.txt"
    elif workspace == 2 or '2':
        root = "/media/nas/nikita/nk_brain_f"
        sc_path = "/media/nas/nikita/nk_brain_f/nk/atlas/hcpmmp1_scale_length.csv"
        lookup_path = "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt"
    elif workspace == 3:
        # HPC subjects
        root = input("Root: ")
        subject_dirs = get_subject_dirs(root)
        for index, directory in enumerate(subject_dirs, start=1):
            print(f"{index}. {directory}")
        choice = int(input("Select a directory by entering its number: "))
        selected_directory = subject_dirs[choice - 1]  # Adjust for zero-based indexing
        sc_path = os.path.join(root, selected_directory, "atlas", "hcpmmp1.csv")
        lookup_path = "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt"
    else:
        raise ValueError("Unknown workspace value provided. Please choose 1, 2, or 3.")

    return root, sc_path, lookup_path


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


def load_node_metrics_as_dataframe(threshold_to_node_csv):
    """
    Given a dictionary mapping {threshold_value: node_csv_path},
    loads all CSVs and concatenates them into one DataFrame with
    columns like:
        ["Label", "Degree Centrality", "Strength",
         "Eigenvector Centrality", "Betweenness Centrality", "Threshold"]

    Returns a single DataFrame containing all thresholds.
    """

    dataframes = []
    for threshold, csv_file in threshold_to_node_csv.items():
        df = pd.read_csv(csv_file)
        # The node CSV has columns:
        #   ["Label", "Degree Centrality", "Strength",
        #    "Eigenvector Centrality", "Betweenness Centrality"]

        # Add a 'Threshold' column so we know which threshold these rows correspond to
        df["Threshold"] = threshold
        dataframes.append(df)

    # Concatenate all per-threshold data into one big DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)
    return combined_df