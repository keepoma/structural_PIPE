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


import os
import pandas as pd


def average_group_file(root_dir, session, file_rel_path, output_file):
    """
    Averages CSV files across subjects for a given session and relative file path.

    Parameters:
      - root_dir: The directory containing subject folders (e.g., "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats")
      - session: The session folder name ("ses_pre" or "ses_post")
      - file_rel_path: The relative path from the subject folder to the CSV file.
                       For example: "connectome_stats/t90/hcpmmp1_minmax_t90_node_metrics.csv"
      - output_file: The full path for the output averaged CSV.

    The function searches for all subject folders starting with "sub-" inside root_dir, reads the CSV
    file for the given session, averages the numeric columns (using the "Label" column for alignment if present),
    and saves the result to output_file.
    """
    # Find subject directories that start with "sub-"
    subject_dirs = [os.path.join(root_dir, d) for d in os.listdir(root_dir)
                    if os.path.isdir(os.path.join(root_dir, d)) and d.startswith("sub-")]

    if not subject_dirs:
        print("No subject directories found in", root_dir)
        return

    dfs = []
    for sub in subject_dirs:
        csv_path = os.path.join(sub, session, file_rel_path)
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path)
                # If there's a Label column, set it as the index so we average by node labels.
                if "Label" in df.columns:
                    df.set_index("Label", inplace=True)
                dfs.append(df)
            except Exception as e:
                print(f"Error reading {csv_path}: {e}")
        else:
            print("File not found:", csv_path)

    if not dfs:
        print("No CSV files found for session", session, "with relative path", file_rel_path)
        return

    # Concatenate dataframes and average by index (node labels) using only numeric columns.
    avg_df = pd.concat(dfs, axis=0).groupby(level=0).mean(numeric_only=True)
    avg_df.reset_index(inplace=True)

    # Save the averaged dataframe
    avg_df.to_csv(output_file, index=False)
    print("Averaged CSV saved to", output_file)


import os
import numpy as np
import pandas as pd

def average_third_column_edges(root_dir, session, file_rel_path, output_file):
    """
    Reads the same edge-metrics CSV from each subject's session folder and averages
    only the 3rd column. The 1st and 2nd columns are kept exactly as in the first
    (reference) CSV, preserving both the column names and row order.

    Parameters:
      root_dir: The parent directory containing sub-XX folders.
      session:  The session folder name (e.g., "ses_pre" or "ses_post").
      file_rel_path: The relative path (from each sub-XX/session folder) to the edges CSV.
      output_file: The full path to the output CSV with averaged 3rd column.
    """

    # Identify subject directories (e.g. sub-11, sub-12, etc.)
    subject_dirs = [
        os.path.join(root_dir, d)
        for d in os.listdir(root_dir)
        if d.startswith("sub-") and os.path.isdir(os.path.join(root_dir, d))
    ]
    if not subject_dirs:
        print(f"No subject directories found in {root_dir}")
        return

    reference_df = None   # Will store the first CSV read (the "reference")
    collected_third_cols = []  # List of arrays for the 3rd column from each subject

    for sub_dir in subject_dirs:
        csv_path = os.path.join(sub_dir, session, file_rel_path)
        if not os.path.exists(csv_path):
            print(f"File not found: {csv_path}")
            continue

        try:
            df = pd.read_csv(csv_path)  # assumes the CSV has a header row
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")
            continue

        # The CSV must have at least 3 columns.
        if df.shape[1] < 3:
            print(f"{csv_path} does not have at least 3 columns. Skipping.")
            continue

        if reference_df is None:
            # Use this first file as our "reference"
            reference_df = df.copy()
            collected_third_cols.append(df.iloc[:, 2].values)
        else:
            # Check if columns match the reference (both names and count)
            if not df.columns.equals(reference_df.columns):
                print(f"Column mismatch in {csv_path}. Skipping.")
                continue
            # Check if row counts match
            if df.shape[0] != reference_df.shape[0]:
                print(f"Row count mismatch in {csv_path}. Skipping.")
                continue
            # Check if first two columns are identical (exact same pairs in same order)
            if not df.iloc[:, :2].equals(reference_df.iloc[:, :2]):
                print(f"First two columns mismatch in {csv_path}. Skipping.")
                continue

            # If everything matches, collect the 3rd column
            collected_third_cols.append(df.iloc[:, 2].values)

    if reference_df is None or not collected_third_cols:
        print("No valid CSVs found to average.")
        return

    # Stack all the 3rd-column arrays and compute the mean row-wise
    stacked = np.vstack(collected_third_cols)  # shape: (num_subjects, num_rows)
    avg_third_col = stacked.mean(axis=0)       # shape: (num_rows,)

    # Place the averaged values back into the 3rd column of the reference DataFrame
    reference_df.iloc[:, 2] = avg_third_col

    # Write out the final CSV (keeps original column names, row order, 1st & 2nd columns)
    reference_df.to_csv(output_file, index=False)
    print("Averaged edges saved to:", output_file)


root_dir = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
output_dir = os.path.join(root_dir, "group_analysis", "t90")
os.makedirs(output_dir, exist_ok=True)

# For ses_pre node metrics:
average_group_file(
    root_dir=root_dir,
    session="ses_pre",
    file_rel_path="connectome_stats/unthresholded/hcpmmp1_minmax_unmodified_node_metrics.csv",
    output_file=os.path.join(output_dir, "unthresholded_node_metrics.csv")
)

# For ses_pre edge metrics:
average_third_column_edges(
    root_dir=root_dir,
    session="ses_pre",
    file_rel_path="connectome_stats/unthresholded/hcpmmp1_minmax_unmodified_top_50_edges.csv",
    output_file=os.path.join(output_dir, "edgesavgNT.csv")
)

# For ses_post node metrics:
average_group_file(
    root_dir=root_dir,
    session="ses_post",
    file_rel_path="connectome_stats/t90/hcpmmp1_invleng_t90_node_metrics.csv",
    output_file=os.path.join(output_dir, "nodes_post_invleng.csv")
)

# For ses_post edge metrics:
average_third_column_edges(
    root_dir=root_dir,
    session="ses_post",
    file_rel_path="connectome_stats/t90/hcpmmp1_invleng_t90_top_50_edges.csv",
    output_file=os.path.join(output_dir, "edges_post_invleng.csv")
)

