import numpy as np
import os
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import pandas as pd
from helpers.j_helpers import get_subject_dirs

def visualize_matrix_weights(matrix, bins=50):
    """
    Loads a matrix and plots a histogram of its values.
    """

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


def visualize_matrix_comparison(
    original_matrix,
    matrix_cut_w_mean,
    matrix_cut_w_1std,
    matrix_cut_w_2std,
    title="Comparison of Thresholded Matrices",
    bins=50
):
    """
    Plots histograms of four different matrices (or 1D arrays) in a 2x2 layout:
      1) The original matrix
      2) matrix_cut_w_mean (cut at mean)
      3) matrix_cut_w_1std (cut at 1*std)
      4) matrix_cut_w_2std (cut at 2*std)
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=16)

    # Top-left: Original
    axs[0, 0].hist(original_matrix.ravel(), bins=bins, edgecolor='black')
    axs[0, 0].set_title("Original")
    axs[0, 0].set_xlabel("Weight")
    axs[0, 0].set_ylabel("Frequency")

    # Top-right: Cut >= Mean
    axs[0, 1].hist(matrix_cut_w_mean.ravel(), bins=bins, edgecolor='black')
    axs[0, 1].set_title("Cut >= Mean")
    axs[0, 1].set_xlabel("Weight")
    axs[0, 1].set_ylabel("Frequency")

    # Bottom-left: Cut >= 1 STD
    axs[1, 0].hist(matrix_cut_w_1std.ravel(), bins=bins, edgecolor='black')
    axs[1, 0].set_title("Cut >= 1 STD")
    axs[1, 0].set_xlabel("Weight")
    axs[1, 0].set_ylabel("Frequency")

    # Bottom-right: Cut >= 2 STD
    axs[1, 1].hist(matrix_cut_w_2std.ravel(), bins=bins, edgecolor='black')
    axs[1, 1].set_title("Cut >= 2 STD")
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


def visualize_matrix(matrix, clip):
    if clip:
        upper_bound = np.percentile(matrix, 99)
        matrix = np.clip(matrix, None, upper_bound)
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, cmap='viridis', square=True)
    plt.title("title")
    plt.xlabel('Region Index')
    plt.ylabel('Region Index')
    plt.show()


def visualize_matrix_side_by_side(matrix):

    # Create a clipped version (using the 99th percentile as the upper bound)
    upper_bound = np.percentile(matrix, 99)
    matrix_clipped = np.clip(matrix, None, upper_bound)

    # Create two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Plot the original (unclipped) matrix
    sns.heatmap(matrix, cmap='viridis', square=True, ax=ax1, cbar_kws={'label': 'Connection Weight'})
    ax1.set_title("Unclipped Matrix")
    ax1.set_xlabel("Region Index")
    ax1.set_ylabel("Region Index")

    # Plot the clipped matrix
    sns.heatmap(matrix_clipped, cmap='viridis', square=True, ax=ax2, cbar_kws={'label': 'Connection Weight'})
    ax2.set_title("Clipped Matrix (99th Percentile)")
    ax2.set_xlabel("Region Index")
    ax2.set_ylabel("Region Index")

    # Set an overall title for the figure
    fig.suptitle("Comparison of Clipped and Unclipped Structural Connectivity Matrix", fontsize=16)

    # Adjust layout to accommodate the overall title
    plt.tight_layout(rect=[0, 0, 1, 0.95])
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


def gather_pre_post_node_data(
        root,
        thresh_folder="t90",
        metric="Degree Centrality",
        average_across_subjects=False
):
    """
    Reads the node-level metrics for all subjects from 'thresh_folder' (e.g., 't90')
    for both ses_pre and ses_post, and merges them by subject + node label (i.e., "Subject"
    and "Label"). If 'average_across_subjects' is True, then the function will group
    by 'Label' and average across subjects, resulting in a single row per node label.
    Otherwise, it returns one row per subject + node combination.

    Returns a DataFrame with columns:
      - if average_across_subjects=False:
         ['Subject', 'Label', 'pre_value', 'post_value']
      - if average_across_subjects=True:
         ['Label', 'pre_value', 'post_value']

    The node metrics filenames are assumed to be something like:
      hcpmmp1_minmax_t90_node_metrics.csv
    in each subject's directory:
      sub-XX/ses_pre/connectome_stats/t90/
      sub-XX/ses_post/connectome_stats/t90/

    Parameters:
      root (str): The root directory containing subject folders.
      thresh_folder (str): The name of the threshold subfolder (e.g., "t90", "t80", or "unthresholded").
      metric (str): The column name of the metric we want to gather (e.g., "Degree Centrality").
      average_across_subjects (bool): If True, the result is grouped by node label, giving one row
                                      per node (averaged across subjects). If False, we keep the
                                      subject-level detail.

    Returns:
      pd.DataFrame or None if no data found.
    """
    pre_records = []
    post_records = []

    # Build the node metrics filename, e.g. hcpmmp1_minmax_t90_node_metrics.csv
    if thresh_folder == "unthresholded":
        node_filename = "hcpmmp1_minmax_unmodified_node_metrics.csv"
    else:
        node_filename = f"hcpmmp1_minmax_{thresh_folder}_node_metrics.csv"

    subjects = get_subject_dirs(root)
    for subj in subjects:
        subj_path = os.path.join(root, subj)

        pre_csv = os.path.join(subj_path, "ses_pre", "connectome_stats", thresh_folder, node_filename)
        post_csv = os.path.join(subj_path, "ses_post", "connectome_stats", thresh_folder, node_filename)

        if not os.path.isfile(pre_csv) or not os.path.isfile(post_csv):
            # Skip if either file is missing
            continue

        df_pre = pd.read_csv(pre_csv)
        df_post = pd.read_csv(post_csv)

        # Keep only the columns we need: Label + metric
        # We'll add Subject so we can later keep or group by it
        df_pre = df_pre[["Label", metric]].copy()
        df_pre["Subject"] = subj

        df_post = df_post[["Label", metric]].copy()
        df_post["Subject"] = subj

        # We'll label them "pre" or "post" if we need them for clarity
        df_pre["Session"] = "pre"
        df_post["Session"] = "post"

        pre_records.append(df_pre)
        post_records.append(df_post)

    if not pre_records or not post_records:
        print("No subjects with both pre and post files found.")
        return None

    # Combine all subjects' pre data, all subjects' post data
    df_pre_all = pd.concat(pre_records, ignore_index=True)
    df_post_all = pd.concat(post_records, ignore_index=True)

    # Rename the metric column to pre_value or post_value
    df_pre_all.rename(columns={metric: "pre_value"}, inplace=True)
    df_post_all.rename(columns={metric: "post_value"}, inplace=True)

    # Merge on Subject + Label to get pre_value and post_value side by side
    merged = pd.merge(
        df_pre_all[["Subject", "Label", "pre_value"]],
        df_post_all[["Subject", "Label", "post_value"]],
        on=["Subject", "Label"],
        how="inner"
    )

    if merged.empty:
        print("No common (Subject, Label) pairs found.")
        return None

    if average_across_subjects:
        # Group by Label and average across subjects
        merged_grouped = merged.groupby("Label", as_index=False).agg({
            "pre_value": "mean",
            "post_value": "mean"
        })
        return merged_grouped  # columns: Label, pre_value, post_value
    else:
        # Return the full detail: columns [Subject, Label, pre_value, post_value]
        return merged


def plot_pre_post_node_lines(merged_df, metric="Degree Centrality", top_n=None):
    """
    Given a DataFrame with columns [Label, pre_value, post_value],
    plot a 'before-after' line plot showing how each node's mean metric changes
    from pre (x=0) to post (x=1).

    Parameters:
      merged_df (pd.DataFrame): Must have columns ['Label', 'pre_value', 'post_value'].
      metric (str): The name of the metric being plotted (used for the y-axis label and title).
      top_n (int or None): If None (default), plot all nodes. Otherwise, plot only the
                           top 'top_n' nodes with the largest absolute difference
                           (|post_value - pre_value|).
    """
    # Compute the difference (post - pre) and absolute difference
    merged_df["diff"] = merged_df["post_value"] - merged_df["pre_value"]
    merged_df["abs_diff"] = merged_df["diff"].abs()

    # If user requested top_n nodes by absolute difference
    if top_n is not None:
        merged_df = merged_df.nlargest(top_n, "abs_diff")

    # Sort by actual diff so lines are drawn from smallest to largest difference
    merged_df.sort_values("diff", inplace=True)

    plt.figure(figsize=(8, max(6, 0.2 * len(merged_df))))  # dynamic height

    for i, row in merged_df.iterrows():
        label = row["Label"]
        pre_val = row["pre_value"]
        post_val = row["post_value"]

        # x=0 for pre, x=1 for post
        plt.plot([0, 1], [pre_val, post_val], marker='o')

        # Label near the post point
        plt.text(1.01, post_val, label, va='center', fontsize=8)

    plt.xlim(-0.1, 1.2)
    plt.xticks([0, 1], ["Pre", "Post"])
    plt.xlabel("Session")
    plt.ylabel(metric)

    n_nodes = len(merged_df)
    title_str = f"{metric} (Mean Across Subjects) Pre vs Post - {n_nodes} nodes"
    if top_n is not None:
        title_str += f" (Top {top_n} by {metric})"
    plt.title(title_str)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    # We'll focus on t90 threshold
    thresh_folder = "t90"
    metric = "Eigenvector Centrality"  # or "Strength", "Degree Centrality",
    # "Eigenvector Centrality", "Betweenness Centrality" etc.

    merged = gather_pre_post_node_data(
        root,
        thresh_folder=thresh_folder,
        metric=metric,
        average_across_subjects=True
    )
    if merged is not None:
        plot_pre_post_node_lines(merged, metric=metric, top_n=5)