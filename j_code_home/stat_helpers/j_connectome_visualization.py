import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import seaborn as sns
import pandas as pd
from helpers.j_helpers import get_subject_dirs
from j_statistical_helpers import shorten_label, extract_sub_ses


node_type_dict = {
    # Visual Cortex
    "L_V1": "Visual Cortex",
    "R_V1": "Visual Cortex",
    "L_V2": "Visual Cortex",
    "R_V2": "Visual Cortex",
    "L_V3": "Visual Cortex",
    "R_V3": "Visual Cortex",

    # frontal motor Cortex
    "R_4": "Frontal Lobe",

    # Somatosensory Cortex
    "L_3a": "Parietal Lobe",

    # Diencephalon
    "L_Thalamus": "Diencephalon",
    "R_Thalamus": "Diencephalon",
    "R_VentralDC": "Diencephalon",

    # Basal Ganglia
    "L_Putamen": "Basal Ganglia",
    "R_Putamen": "Basal Ganglia",
    "L_Caudate": "Basal Ganglia",
    "R_Caudate": "Basal Ganglia",

    # Brainstem
    "Brain-Stem": "Brainstem",

    # Cerebellum
    "L_Cerebellum": "Cerebellum",
    "R_Cerebellum": "Cerebellum"
}

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


def visualize_top_n_nodes_and_their_edges(
    node_metrics_csv,
    edge_csv,
    node_metric='Strength',
    n=10,
    layout='circular',
    node_type_dict=None,
    offset=0.07,
    save_figure=False
):
    """
    1. Reads a node metrics CSV containing columns like:
         [Label, Degree Centrality, Strength, Eigenvector Centrality, Betweenness Centrality].
    2. Reads an edges CSV containing columns:
         [Source, Target, Weight].
    3. Selects the top N nodes by the specified node_metric.
    4. Filters edges to only include those connecting the top N nodes among themselves.
    5. Builds a subgraph from these nodes and edges and visualizes it, with node color
       determined by their anatomical category (node_type_dict) and labels offset above nodes.

    Parameters:
      node_metrics_csv (str): Path to the CSV with node-level metrics.
      edge_csv (str): Path to the CSV with edges [Source, Target, Weight].
      node_metric (str): The node-level metric by which to rank (default='Strength').
      n (int): How many top nodes to include (default=10).
      layout (str): Which NetworkX layout to use ('circular', 'kamada_kawai', 'spring', etc.).
      node_type_dict (dict): Mapping {node_label: category}, e.g. 'L_V1' -> 'Visual Cortex'.
      offset (float): Vertical offset for node labels to reduce overlap (default=0.07).
    """

    # Example category-to-color mapping for your chosen categories
    category_color_map = {
        "Visual Cortex": "lightblue",
        "Motor Cortex": "red",
        "Somatosensory Cortex": "olive",
        "Basal Ganglia": "lightgreen",
        "Diencephalon": "pink",
        "Brainstem": "orange",
        "Cerebellum": "lightgray"
    }

    # 1. Read the node metrics
    node_df = pd.read_csv(node_metrics_csv)
    node_df = node_df.sort_values(node_metric, ascending=False)
    top_n_nodes = node_df.head(n)['Label'].tolist()

    # 2. Read the edges
    edge_df = pd.read_csv(edge_csv)

    # 3. Filter edges
    sub_edges = edge_df[
        edge_df['Source'].isin(top_n_nodes) & edge_df['Target'].isin(top_n_nodes)
    ]

    # 4. Build the subgraph
    G = nx.Graph()

    if node_type_dict is None:
        # If no dictionary is provided, default to "Visual Cortex" or any fallback
        node_type_dict = {}

    # Add top N nodes
    for _, row in node_df.head(n).iterrows():
        label = row['Label']
        category = node_type_dict.get(label, "Visual Cortex")  # fallback if not in dict
        G.add_node(label, **row.to_dict(), category=category)

    # Add edges
    for _, row in sub_edges.iterrows():
        source = row['Source']
        target = row['Target']
        weight = float(row['Weight'])
        G.add_edge(source, target, weight=weight)

    # 5. Choose layout
    if layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G, weight='weight')
    elif layout == 'spring':
        pos = nx.spring_layout(G, seed=42)
    else:
        pos = nx.spring_layout(G, seed=42)

    plt.figure(figsize=(7, 7))

    # Scale node sizes by the chosen metric
    metric_values = node_df[node_metric]
    max_val = metric_values.max() if len(metric_values) > 0 else 1
    node_sizes = []
    for node in G.nodes():
        node_val = G.nodes[node].get(node_metric, 1)
        node_sizes.append(300 + 2500 * (node_val / max_val))

    # Determine node colors by category
    node_colors = []
    for node in G.nodes():
        cat = G.nodes[node].get('category', 'Visual Cortex')
        color = category_color_map.get(cat, 'lightblue')  # fallback color
        node_colors.append(color)

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors='black'
    )

    # Scale edge widths
    edges = G.edges(data=True)
    if edges:
        weights = [edge[2]['weight'] for edge in edges]
        w_min, w_max = min(weights), max(weights)
        if w_min == w_max:
            scaled_widths = [3]*len(weights)
        else:
            scaled_widths = [
                1 + (5 - 1) * ((w - w_min) / (w_max - w_min))
                for w in weights
            ]
    else:
        scaled_widths = []

    nx.draw_networkx_edges(
        G, pos,
        width=scaled_widths,
        edge_color='gray'
    )

    # Offset labels above nodes
    label_positions = {}
    display_labels = {}
    for node, (x, y) in pos.items():
        label_positions[node] = (x, y + offset)
        display_labels[node] = shorten_label(node)
    nx.draw_networkx_labels(
        G,
        label_positions,
        labels=display_labels,
        font_size=9
    )
    # Build legend handles from category_color_map
    legend_handles = []
    for cat, color in category_color_map.items():
        legend_handles.append(
            Line2D([0], [0], marker='o', color='w', label=cat,
                   markerfacecolor=color, markersize=10, markeredgecolor='black')
        )

    plt.legend(
        handles=legend_handles,
        loc='upper right',
        bbox_to_anchor=(1.1, 0.98),  # shift legend slightly outside the axes to the right
        borderaxespad=0.,
        fontsize=9
    )

    subj_ses = extract_sub_ses(node_metrics_csv)

    title_text = f"Top {n} Nodes by {node_metric} with Interconnecting Edges for {subj_ses}"
    plt.title(title_text)
    plt.axis('off')

    if save_figure:
        save_dir = "/Users/nikitakaruzin/Desktop/Research/Picht/Lab Report/figs"
        os.makedirs(save_dir, exist_ok=True)  # ensure directory exists

        # Create a safe filename from the title by replacing spaces with underscores
        safe_filename = title_text.replace(' ', '_') + ".png"
        save_path = os.path.join(save_dir, safe_filename)

        # Save with higher resolution (dpi=300) and tight bounding box
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()


def visualize_two_files_side_by_side(
    node_metrics_csv_pre,
    edge_csv_pre,
    node_metrics_csv_post,
    edge_csv_post,
    node_metric='Strength',
    n=10,
    layout='circular',
    node_type_dict=None,
    offset=0.07,
    save_figure=False
):
    """
    Creates a side-by-side figure with two subplots:
      - Left: the top N nodes/edges from (node_metrics_csv_pre, edge_csv_pre)
      - Right: the top N nodes/edges from (node_metrics_csv_post, edge_csv_post)
    This allows a direct visual comparison (e.g., pre- vs post-therapy).

    Parameters:
      node_metrics_csv_pre (str): CSV for node-level metrics (pre).
      edge_csv_pre (str): CSV for edges [Source, Target, Weight] (pre).
      node_metrics_csv_post (str): CSV for node-level metrics (post).
      edge_csv_post (str): CSV for edges (post).
      node_metric (str): Which metric to rank by (default='Strength').
      n (int): Number of top nodes to display (default=10).
      layout (str): Layout name ('circular', 'kamada_kawai', 'spring', etc.).
      node_type_dict (dict): Mapping for color-coding categories (optional).
      offset (float): Vertical offset for node labels (default=0.07).
      save_figure (bool): If True, saves the entire side-by-side figure.

    Returns:
      None
    """

    # Common color map
    category_color_map = {
        "Visual Cortex": "lightblue",
        "Frontal Lobe": "gold",
        "Parietal Lobe": "olive",
        "Basal Ganglia": "lightgreen",
        "Diencephalon": "pink",
        "Brainstem": "orange",
        "Cerebellum": "lightgray"
    }

    # --------------------
    # 1) Build PRE Graph
    # --------------------
    df_pre_nodes = pd.read_csv(node_metrics_csv_pre).sort_values(node_metric, ascending=False)
    top_n_pre = df_pre_nodes.head(n)['Label'].tolist()

    df_pre_edges = pd.read_csv(edge_csv_pre)
    sub_pre_edges = df_pre_edges[
        df_pre_edges['Source'].isin(top_n_pre) & df_pre_edges['Target'].isin(top_n_pre)
    ]

    G_pre = nx.Graph()
    if node_type_dict is None:
        node_type_dict = {}

    for _, row in df_pre_nodes.head(n).iterrows():
        label = row['Label']
        category = node_type_dict.get(label, "Visual Cortex")
        G_pre.add_node(label, **row.to_dict(), category=category)

    for _, row in sub_pre_edges.iterrows():
        G_pre.add_edge(row['Source'], row['Target'], weight=float(row['Weight']))

    if layout == 'circular':
        pos_pre = nx.circular_layout(G_pre)
    elif layout == 'kamada_kawai':
        pos_pre = nx.kamada_kawai_layout(G_pre, weight='weight')
    elif layout == 'spring':
        pos_pre = nx.spring_layout(G_pre, seed=42)
    else:
        pos_pre = nx.spring_layout(G_pre, seed=42)

    # --------------------
    # 2) Build POST Graph
    # --------------------
    df_post_nodes = pd.read_csv(node_metrics_csv_post).sort_values(node_metric, ascending=False)
    top_n_post = df_post_nodes.head(n)['Label'].tolist()

    df_post_edges = pd.read_csv(edge_csv_post)
    sub_post_edges = df_post_edges[
        df_post_edges['Source'].isin(top_n_post) & df_post_edges['Target'].isin(top_n_post)
    ]

    G_post = nx.Graph()
    for _, row in df_post_nodes.head(n).iterrows():
        label = row['Label']
        category = node_type_dict.get(label, "Visual Cortex")
        G_post.add_node(label, **row.to_dict(), category=category)

    for _, row in sub_post_edges.iterrows():
        G_post.add_edge(row['Source'], row['Target'], weight=float(row['Weight']))

    if layout == 'circular':
        pos_post = nx.circular_layout(G_post)
    elif layout == 'kamada_kawai':
        pos_post = nx.kamada_kawai_layout(G_post, weight='weight')
    elif layout == 'spring':
        pos_post = nx.spring_layout(G_post, seed=42)
    else:
        pos_post = nx.spring_layout(G_post, seed=42)

    # --------------------
    # 3) Create figure with 2 subplots
    # --------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle(f"Top {n} Nodes by averaged {node_metric} in five pre-TMS Connectomes", fontsize=16)

    # Helper function to draw one subgraph on a given axis
    def draw_subgraph(ax, G, pos, df_nodes):
        # Node sizes
        metric_values = df_nodes[node_metric]
        max_val = metric_values.max() if len(metric_values) > 0 else 1
        node_sizes = []
        for node in G.nodes():
            node_val = G.nodes[node].get(node_metric, 1)
            node_sizes.append(300 + 2500 * (node_val / max_val))

        # Node colors
        node_colors = []
        for node in G.nodes():
            cat = G.nodes[node].get('category', 'Visual Cortex')
            node_colors.append(category_color_map.get(cat, 'lightblue'))

        # Edges
        edges_data = G.edges(data=True)
        if edges_data:
            weights = [ed[2]['weight'] for ed in edges_data]
            w_min, w_max = min(weights), max(weights)
            if w_min == w_max:
                scaled_widths = [3]*len(weights)
            else:
                scaled_widths = [
                    1 + (5 - 1)*((w - w_min)/(w_max - w_min))
                    for w in weights
                ]
        else:
            scaled_widths = []

        # Draw
        nx.draw_networkx_nodes(G, pos, ax=ax,
                               node_size=node_sizes,
                               node_color=node_colors,
                               edgecolors='black')
        nx.draw_networkx_edges(G, pos, ax=ax,
                               width=scaled_widths,
                               edge_color='gray')

        # Offset labels
        label_positions = {}
        display_labels = {}
        for node, (x, y) in pos.items():
            label_positions[node] = (x, y + offset)
            display_labels[node] = shorten_label(node)

        nx.draw_networkx_labels(G, label_positions, ax=ax,
                                labels=display_labels, font_size=9)
        ax.axis('off')

    # Draw the "pre" graph on ax1
    draw_subgraph(ax1, G_pre, pos_pre, df_pre_nodes)
    subj_ses_pre = extract_sub_ses(node_metrics_csv_pre)
    ax1.set_title(f"Connectome scaled by inverse length", fontsize=14)

    # Draw the "post" graph on ax2
    draw_subgraph(ax2, G_post, pos_post, df_post_nodes)
    subj_ses_post = extract_sub_ses(node_metrics_csv_post)
    ax2.set_title(f"Unscaled connectome", fontsize=14)

    # Build legend handles once
    legend_handles = []
    for cat, color in category_color_map.items():
        legend_handles.append(
            Line2D([0], [0], marker='o', color='w', label=cat,
                   markerfacecolor=color, markersize=10, markeredgecolor='black')
        )

    # Place a single legend for the figure (outside the subplots)
    fig.legend(
        handles=legend_handles,
        loc='upper right',
        bbox_to_anchor=(0.99, 0.99),
        borderaxespad=0.,
        fontsize=9
    )

    # 4) Optionally save figure
    if save_figure:
        save_dir = "/Users/nikitakaruzin/Desktop/Research/Picht/Lab Report/figs"
        os.makedirs(save_dir, exist_ok=True)
        # Create a filename from the suptitle
        safe_filename = f"Side_by_side_Top_{n}_{node_metric}_invlength.png".replace(' ', '_')
        save_path = os.path.join(save_dir, safe_filename)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Side-by-side figure saved to {save_path}")

    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-1.2, 1.2)
    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(-1.2, 1.2)
    plt.tight_layout()

    plt.show()



if __name__ == "__main__":
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    # We'll focus on t90 threshold
    thresh_folder = "t90"
    metric = "Strength"  # or "Strength", "Degree Centrality", Eigenvector Centrality
    nodes_pre_minmax = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/nodes_pre_minmax.csv")
    edges_pre_minmax = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/edges_pre_minmax.csv")

    nodes_post_minmax = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/nodes_post_minmax.csv")
    edges_post_minmax = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/edges_post_minmax.csv")

    visualize_two_files_side_by_side(
        nodes_pre_minmax,
        edges_pre_minmax,
        nodes_post_minmax,
        edges_post_minmax,
        node_metric='Strength',
        n=15,
        layout='circular',
        node_type_dict=node_type_dict,
        offset=0,
        save_figure=False
    )

    nodes_pre_invlenght = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/nodes_pre_invleng.csv")
    edges_pre_invlenght = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/edges_pre_invleng.csv")

    nodes_post_invlenght = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/nodes_post_invleng.csv")
    edges_post_invlenght = ("/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/edges_post_invleng.csv")
    visualize_two_files_side_by_side(
        nodes_pre_invlenght,
        edges_pre_invlenght,
        nodes_pre_minmax,
        edges_pre_minmax,
        node_metric='Strength',
        n=15,
        layout='circular',
        node_type_dict=node_type_dict,
        offset=0,
        save_figure=False
    )


    """
    merged = gather_pre_post_node_data(
        root,
        thresh_folder=thresh_folder,
        metric=metric,
        average_across_subjects=True
    )  
    if merged is not None:
        plot_pre_post_node_lines(merged, metric=metric, top_n=5)
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    # Load CSVs
    df_unthresh = pd.read_csv(
        '/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/unthresholded/unthresholded_node_metrics.csv'
    )
    df_thresh = pd.read_csv(
        '/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/group_analysis/t90/nodes_pre_minmax.csv'
    )

    # Index DataFrames by "Label"
    node_col = "Label"
    df_unthresh_lookup = df_unthresh.set_index(node_col)
    df_thresh_lookup = df_thresh.set_index(node_col)

    # Define the metrics and setup subplots
    metrics = [
        "Degree Centrality",
        "Strength",
        "Eigenvector Centrality",
        "Betweenness Centrality"
    ]

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    axes = axes.flatten()

    # Choose two shades of blue for plotting
    color_unth = "lightblue"
    color_thresh = "cornflowerblue"

    # Loop over each metric
    for i, metric in enumerate(metrics):
        ax = axes[i]

        # Sort the unthresholded DataFrame in descending order by the current metric
        df_sorted = df_unthresh.sort_values(by=metric, ascending=False)
        all_nodes = df_sorted[node_col].values

        # Pick 20 nodes at regular intervals across the sorted data
        indices = np.linspace(0, len(all_nodes) - 1, 20, dtype=int)
        selected_nodes = all_nodes[indices]

        # Get unthresholded and thresholded values for these nodes
        unthresh_values = [
            df_unthresh_lookup.at[node, metric] if node in df_unthresh_lookup.index else 0
            for node in selected_nodes
        ]
        thresh_values = [
            df_thresh_lookup.at[node, metric] if node in df_thresh_lookup.index else 0
            for node in selected_nodes
        ]

        # Plot side-by-side horizontal bars for each selected node
        y_pos = np.arange(len(selected_nodes))
        bar_height = 0.4

        ax.barh(
            y_pos - bar_height / 2,
            unthresh_values,
            height=bar_height,
            color=color_unth,
            label="Unthresholded"
        )
        ax.barh(
            y_pos + bar_height / 2,
            thresh_values,
            height=bar_height,
            color=color_thresh,
            label="Thresholded"
        )

        # Subplot styling
        ax.set_yticks(y_pos)
        ax.set_yticklabels(selected_nodes, fontsize=7)
        ax.invert_yaxis()  # Largest at the top
        ax.set_xlabel("Value", fontsize=8)
        ax.set_title(f"Sampled 20 Nodes by {metric}", fontsize=9)

        # Show legend on the first subplot (or on each if desired)
        if i == 0:
            ax.legend(fontsize=7)

        # Apply custom tick formatting for the "Strength" metric
        if metric == "Strength":
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x / 1000:.0f}k' if x >= 1000 else f'{x:.0f}'))

    # Main title and layout adjustments
    plt.suptitle(
        "Comparing Unthresholded vs. Thresholded (T90)\nSampling 20 Nodes at Regular Intervals per Metric",
        fontsize=10
    )
    plt.tight_layout()
    plt.show()
