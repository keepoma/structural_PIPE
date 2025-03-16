import os
import numpy as np
import networkx as nx
import pandas as pd
import re


def get_subject_dirs(root):
    """
    Return a sorted list of subject directories excluding a specific folder
    """

    subject_dirs = sorted([
        os.path.join(root, d)
        for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)) and d not in ("group_analysis", "logs")
    ])
    return subject_dirs



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


def shorten_label(label):
    """
    If 'label' contains an underscore, we split on the first underscore.
    If the suffix (the part after the underscore) is longer than 3 characters,
    we truncate it to the first 3. Otherwise, we leave it as is.

    Examples:
      - "L_Cerebellum" -> "L_Cer"
      - "L_V2" -> "L_V2" (since 'V2' is only 2 chars)
      - "Brainstem" -> "Brainstem" (no underscore)
    """
    if "_" in label:
        prefix, suffix = label.split("_", 1)
        if len(suffix) > 3:
            suffix = suffix[:3]
        return prefix + "_" + suffix
    else:
        return label


def extract_sub_ses(path):
    """
    Given a path like:
      /Users/.../sub-11/ses_post/connectome_stats/...
    This function returns:
      'Subject 11 Post-Therapy'
    If it's sub-11/ses_pre, it returns:
      'Subject 11 Pre-Therapy'
    """
    parts = path.split(os.sep)

    subject_part = None
    session_part = None

    # Identify the parts that start with "sub-" or "ses_"
    for part in parts:
        if part.startswith("sub-"):
            subject_part = part
        elif part.startswith("ses_"):
            session_part = part

    if not subject_part or not session_part:
        return None  # Could not find both parts

    # Extract the subject number from sub-XX (e.g., "sub-11" -> "11")
    # We'll use a regex to handle numeric IDs. If your IDs have letters, adjust accordingly.
    match = re.match(r"sub-(\d+)", subject_part)
    if match:
        subject_num = match.group(1)
    else:
        # Fallback: remove "sub-" if the pattern is guaranteed, or leave as-is if uncertain
        subject_num = subject_part.replace("sub-", "")

    # Map session to a descriptive phrase
    if session_part == "ses_pre":
        therapy_str = "Pre-TMS"
    elif session_part == "ses_post":
        therapy_str = "Post-TMS"
    else:
        therapy_str = session_part  # Fallback if you have other sessions

    return f"Subject {subject_num} {therapy_str}"

