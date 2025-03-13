import numpy as np
import os
import networkx as nx
import csv
import pandas as pd
from scipy import stats
from stat_helpers.j_statistical_helpers import (lookup_dictionary, threshold_matrix_by_weight,
                                         threshold_matrix_by_clipping, create_graph,
                                         load_node_metrics_as_dataframe, get_subject_dirs)
from scipy.stats import wilcoxon


def threshold_and_save_matrix_by_top_percent(matrix, top_percent, matrix_path):
    """
    Thresholds the connectivity matrix to keep only the top 'top_percent'
    of nonzero connection values.

    Parameters:
      matrix (np.array): The connectivity matrix.
      top_percent (float): The percentage of strongest values to retain (e.g., 10 or 20).

    Returns:
      thresholded_matrix (np.array): The matrix with values below the threshold set to 0.
      threshold (float): The threshold value used.
    """
    # Extract nonzero values to avoid bias from zeros
    nonzero_values = matrix[matrix > 0]
    if nonzero_values.size == 0:
        threshold = 0
    else:
        # E.g: For top 10%, compute the 90th percentile (100 - 10 = 90)
        threshold = np.percentile(nonzero_values, 100 - top_percent)

    # Threshold the matrix: keep values above or equal to the threshold
    thresholded_matrix = np.where(matrix >= threshold, matrix, 0)

    # Construct new file name: add suffix before extension
    base, ext = os.path.splitext(matrix_path)
    new_file_path = f"{base}_top{int(top_percent)}{ext}"

    # Save the thresholded matrix
    np.savetxt(new_file_path, thresholded_matrix, delimiter=',')
    print(f"Saved thresholded matrix (top {top_percent}%) to: {new_file_path}")
    print(f"Threshold value used: {threshold}")

    return thresholded_matrix


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


def compute_metrics_for_prethresholded_matrix(matrix, lookup_path, output_dir, base_filename, binarize=False,
                                              overwrite=False):
    """
    Computes connectivity metrics from a weighted connectivity matrix that is already thresholded.

    Parameters:
      matrix (np.array): The already thresholded weighted connectivity matrix.
      lookup_path (str): Path to the lookup file for node labels.
      output_dir (str): Directory to save the CSV files.
      base_filename (str): Base filename (without extension) of the source matrix.
      binarize (bool): If True, binarizes the matrix after thresholding.
      overwrite (bool): If True, overwrites existing CSV files.

    Saves:
      - Node-level metrics CSV: base_filename + '_prethresh_node_metrics.csv'
      - Global metrics CSV: base_filename + '_prethresh_global_metrics.csv'

    Returns:
      tuple: (node_metrics_csv_path, global_metrics_csv_path)
    """

    # Optionally binarize the matrix
    if binarize:
        matrix = np.where(matrix > 0, 1.0, 0.0)

    # Create the graph from the matrix
    G = create_graph(matrix)

    # Compute node-level metrics.
    degree_centrality = nx.degree_centrality(G)
    strength = {
        node: sum(data.get('weight', 1) for _, _, data in G.edges(node, data=True))
        for node in G.nodes()
    }
    eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
    betweenness_centrality = nx.betweenness_centrality(G, weight='weight')

    # Compute global metrics.
    grc = global_reaching_centrality(G, centrality_func=nx.degree_centrality)
    global_efficiency = nx.global_efficiency(G)
    local_efficiency = nx.local_efficiency(G)

    # Prepare file names.
    node_csv_file = os.path.join(output_dir, f"{base_filename}_node_metrics.csv")
    global_csv_file = os.path.join(output_dir, f"{base_filename}_global_metrics.csv")

    # Check for file existence.
    if (os.path.isfile(node_csv_file) and os.path.isfile(global_csv_file)) and not overwrite:
        print("CSV files already exist; skipping computation.")
        return node_csv_file, global_csv_file

    # Load lookup dictionary for node labels.
    lookup = lookup_dictionary(lookup_path)

    # Save node-level metrics to CSV.
    with open(node_csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Label", "Degree Centrality", "Strength", "Eigenvector Centrality", "Betweenness Centrality"])
        for node in G.nodes():
            label = lookup.get(node, f"Node_{node}")
            d = degree_centrality.get(node, 0)
            s = strength.get(node, 0)
            e = eigenvector_centrality.get(node, 0)
            b = betweenness_centrality.get(node, 0)
            writer.writerow([label, f"{d:.4f}", f"{s:.4f}", f"{e:.4f}", f"{b:.4f}"])

    # Save global metrics to CSV.
    with open(global_csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Global Reaching Centrality", "Global Efficiency", "Local Efficiency"])
        writer.writerow([f"{grc:.4f}", f"{global_efficiency:.4f}", f"{local_efficiency:.4f}"])

    print(f"Saved node-level metrics to: {node_csv_file}")
    print(f"Saved global metrics to: {global_csv_file}")

    return node_csv_file, global_csv_file


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


def get_metric_filenames(thresh_folder):
    """
    Returns the node metrics and global metrics filenames for a given threshold folder.
    For example:
      - unthresholded -> hcpmmp1_minmax_unmodified_node_metrics.csv, hcpmmp1_minmax_unmodified_global_metrics.csv
      - t90          -> hcpmmp1_minmax_t90_node_metrics.csv, hcpmmp1_minmax_t90_global_metrics.csv
      - t80          -> hcpmmp1_minmax_t80_node_metrics.csv, hcpmmp1_minmax_t80_global_metrics.csv
    """
    if thresh_folder == "unthresholded":
        node_file = "hcpmmp1_minmax_unmodified_node_metrics.csv"
        global_file = "hcpmmp1_minmax_unmodified_global_metrics.csv"
    else:
        # For t80, t90, etc.
        node_file = f"hcpmmp1_minmax_{thresh_folder}_node_metrics.csv"
        global_file = f"hcpmmp1_minmax_{thresh_folder}_global_metrics.csv"
    return node_file, global_file


def run_averaged_nodes_statistical_tests(root, thresh_folder):
    """
    Loads node/global metrics from the specified threshold folder (unthresholded/t90/t80),
    averages node-level metrics across nodes, collects global metrics, and runs Wilcoxon tests.

    Saves results in 'group_analysis/<thresh_folder>/'.
    """
    group_analysis_dir = os.path.join(root, "group_analysis", thresh_folder)
    os.makedirs(group_analysis_dir, exist_ok=True)

    subjects = get_subject_dirs(root)
    subjects.sort()

    node_metrics_pre = {}
    node_metrics_post = {}
    global_metrics_pre = {}
    global_metrics_post = {}

    # Build the correct file names for node/global metrics
    node_filename, global_filename = get_metric_filenames(thresh_folder)

    for subj in subjects:
        subj_path = os.path.join(root, subj)
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_path, session, "connectome_stats", thresh_folder)
            node_file = os.path.join(session_dir, node_filename)
            global_file = os.path.join(session_dir, global_filename)

            try:
                df_node = pd.read_csv(node_file)
            except Exception as e:
                print(f"Error loading {node_file}: {e}")
                continue

            # Average across nodes
            avg_degree = df_node["Degree Centrality"].mean()
            avg_strength = df_node["Strength"].mean()
            avg_eigen = df_node["Eigenvector Centrality"].mean()
            avg_betweenness = df_node["Betweenness Centrality"].mean()

            try:
                df_global = pd.read_csv(global_file)
            except Exception as e:
                print(f"Error loading {global_file}: {e}")
                continue

            global_reaching = df_global["Global Reaching Centrality"].iloc[0]
            global_eff = df_global["Global Efficiency"].iloc[0]
            local_eff = df_global["Local Efficiency"].iloc[0]

            if session == 'ses_pre':
                node_metrics_pre[subj] = {
                    "Degree Centrality": avg_degree,
                    "Strength": avg_strength,
                    "Eigenvector Centrality": avg_eigen,
                    "Betweenness Centrality": avg_betweenness
                }
                global_metrics_pre[subj] = {
                    "Global Reaching Centrality": global_reaching,
                    "Global Efficiency": global_eff,
                    "Local Efficiency": local_eff
                }
            else:
                node_metrics_post[subj] = {
                    "Degree Centrality": avg_degree,
                    "Strength": avg_strength,
                    "Eigenvector Centrality": avg_eigen,
                    "Betweenness Centrality": avg_betweenness
                }
                global_metrics_post[subj] = {
                    "Global Reaching Centrality": global_reaching,
                    "Global Efficiency": global_eff,
                    "Local Efficiency": local_eff
                }

    common_subjects = sorted(set(node_metrics_pre.keys()).intersection(node_metrics_post.keys()))
    if not common_subjects:
        print("No common subjects with both pre and post metrics found.")
        return None, None

    node_metrics_list = ["Degree Centrality", "Strength", "Eigenvector Centrality", "Betweenness Centrality"]
    results_node = {}
    for metric in node_metrics_list:
        pre_vals = [node_metrics_pre[subj][metric] for subj in common_subjects]
        post_vals = [node_metrics_post[subj][metric] for subj in common_subjects]
        try:
            stat, p = wilcoxon(pre_vals, post_vals)
            results_node[metric] = (stat, p)
        except Exception as e:
            results_node[metric] = f"Error: {e}"

    global_metrics_list = ["Global Reaching Centrality", "Global Efficiency", "Local Efficiency"]
    results_global = {}
    for metric in global_metrics_list:
        pre_vals = [global_metrics_pre[subj][metric] for subj in common_subjects]
        post_vals = [global_metrics_post[subj][metric] for subj in common_subjects]
        try:
            stat, p = wilcoxon(pre_vals, post_vals)
            results_global[metric] = (stat, p)
        except Exception as e:
            results_global[metric] = f"Error: {e}"

    # Save results
    node_results_df = pd.DataFrame(
        [{"Metric": metric, "Statistic": stat, "P_value": p} for metric, (stat, p) in results_node.items()]
    )
    node_outfile = os.path.join(group_analysis_dir, f"wilcoxon_averaged_node_metrics_results_{thresh_folder}.csv")
    node_results_df.to_csv(node_outfile, index=False)
    print(f"Averaged node-level Wilcoxon test results saved to: {node_outfile}")

    global_results_df = pd.DataFrame(
        [{"Metric": metric, "Statistic": stat, "P_value": p} for metric, (stat, p) in results_global.items()]
    )
    global_outfile = os.path.join(group_analysis_dir, f"wilcoxon_averaged_global_metrics_results_{thresh_folder}.csv")
    global_results_df.to_csv(global_outfile, index=False)
    print(f"Averaged global metrics Wilcoxon test results saved to: {global_outfile}")

    return results_node, results_global


def run_node_level_wilcoxon(root, thresh_folder):
    """
    For each subject, loads the node metrics CSV from the specified threshold folder
    (t80, t90, or unthresholded) for both pre/post sessions, merges on 'Label', and
    runs Wilcoxon tests for each metric per node. Saves results in group_analysis/<thresh_folder>.
    """
    group_analysis_dir = os.path.join(root, "group_analysis", thresh_folder)
    os.makedirs(group_analysis_dir, exist_ok=True)

    subjects = get_subject_dirs(root)
    subjects.sort()

    data = {}
    node_filename, _ = get_metric_filenames(thresh_folder)  # Only need the node filename.

    for subj in subjects:
        subj_path = os.path.join(root, subj)
        pre_csv = os.path.join(subj_path, "ses_pre", "connectome_stats", thresh_folder, node_filename)
        post_csv = os.path.join(subj_path, "ses_post", "connectome_stats", thresh_folder, node_filename)
        try:
            df_pre = pd.read_csv(pre_csv)
            df_post = pd.read_csv(post_csv)
        except Exception as e:
            print(f"Skipping subject {subj}: {e}")
            continue

        for _, row in df_pre.iterrows():
            label = row["Label"]
            deg_pre = row["Degree Centrality"]
            str_pre = row["Strength"]
            eig_pre = row["Eigenvector Centrality"]
            btw_pre = row["Betweenness Centrality"]

            df_post_node = df_post[df_post["Label"] == label]
            if df_post_node.empty:
                print(f"Subject {subj}: Node {label} not found in post CSV. Skipping.")
                continue
            post_row = df_post_node.iloc[0]
            deg_post = post_row["Degree Centrality"]
            str_post = post_row["Strength"]
            eig_post = post_row["Eigenvector Centrality"]
            btw_post = post_row["Betweenness Centrality"]

            if label not in data:
                data[label] = {
                    "Degree Centrality": {"pre": [], "post": []},
                    "Strength": {"pre": [], "post": []},
                    "Eigenvector Centrality": {"pre": [], "post": []},
                    "Betweenness Centrality": {"pre": [], "post": []},
                }
            data[label]["Degree Centrality"]["pre"].append(deg_pre)
            data[label]["Degree Centrality"]["post"].append(deg_post)
            data[label]["Strength"]["pre"].append(str_pre)
            data[label]["Strength"]["post"].append(str_post)
            data[label]["Eigenvector Centrality"]["pre"].append(eig_pre)
            data[label]["Eigenvector Centrality"]["post"].append(eig_post)
            data[label]["Betweenness Centrality"]["pre"].append(btw_pre)
            data[label]["Betweenness Centrality"]["post"].append(btw_post)

    results = []
    metrics = ["Degree Centrality", "Strength", "Eigenvector Centrality", "Betweenness Centrality"]
    for label, metric_data in data.items():
        for metric in metrics:
            pre_vals = metric_data[metric]["pre"]
            post_vals = metric_data[metric]["post"]
            if len(pre_vals) < 2:
                res = {"Label": label, "Metric": metric, "Statistic": None, "P_value": None,
                       "Comment": "Insufficient data"}
            else:
                try:
                    stat, p = wilcoxon(pre_vals, post_vals)
                    res = {"Label": label, "Metric": metric, "Statistic": stat, "P_value": p, "Comment": ""}
                except Exception as e:
                    res = {"Label": label, "Metric": metric, "Statistic": None, "P_value": None, "Comment": str(e)}
            results.append(res)

    results_df = pd.DataFrame(results)
    output_file = os.path.join(group_analysis_dir, f"wilcoxon_node_level_metrics_results_{thresh_folder}.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Wilcoxon node-level test results saved to: {output_file}")

    return results_df


def check_all_p_values(group_analysis_dir, p_threshold=0.05):
    """
    Loops through all subdirectories in 'group_analysis_dir', reads each CSV file,
    and checks if there is a 'P_value' column. If so, prints all rows with p-value < p_threshold.

    Parameters:
      group_analysis_dir (str): Path to the group_analysis directory.
      p_threshold (float): Significance threshold (default=0.05).
    """
    # List all subdirectories in group_analysis_dir
    subdirs = [d for d in os.listdir(group_analysis_dir)
               if os.path.isdir(os.path.join(group_analysis_dir, d))]
    subdirs.sort()

    for subdir in subdirs:
        subdir_path = os.path.join(group_analysis_dir, subdir)
        print(f"\n=== Checking directory: {subdir_path} ===")

        # Loop over all CSV files in this subdirectory
        for file_name in os.listdir(subdir_path):
            if file_name.endswith(".csv"):
                csv_path = os.path.join(subdir_path, file_name)
                try:
                    df = pd.read_csv(csv_path)
                except Exception as e:
                    print(f"Error reading {csv_path}: {e}")
                    continue

                # Check if 'P_value' column exists
                if "P_value" not in df.columns:
                    print(f"Skipping {csv_path}, no 'P_value' column found.")
                    continue

                # Filter rows with p-value < p_threshold
                sig_rows = df[df["P_value"] < p_threshold]
                if sig_rows.empty:
                    print(f"{file_name}: no p-values < {p_threshold} found.")
                else:
                    print(f"{file_name}: Found {len(sig_rows)} rows with p-value < {p_threshold}:")
                    # Print relevant info from each row
                    for idx, row in sig_rows.iterrows():
                        # We check if columns exist before printing
                        label_str = row["Label"] if "Label" in df.columns else "N/A"
                        metric_str = row["Metric"] if "Metric" in df.columns else "N/A"
                        p_val_str = row["P_value"]
                        print(f"  Row {idx}: Label={label_str}, Metric={metric_str}, p={p_val_str}")


if __name__ == '__main__':
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"

    # List of threshold folders we want to process.
    """
    threshold_folders = ["unthresholded", "t90", "t80"]
    for thresh in threshold_folders:
        print(f"\nProcessing averaged metrics for threshold folder: {thresh}")
        run_averaged_nodes_statistical_tests(root, thresh)
        run_node_level_wilcoxon(root, thresh)
    """
    check_all_p_values("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis", p_threshold=0.05)


