import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from helpers.j_helpers import (get_subject_paths_wo_ses)


def plot_mean_fa_along_nodes(subject_dirs, track_base):
    """
    Aggregates node-wise mean FA values across subjects for the left and right
    resampled tracts and creates a line plot with error bars.

    The x-axis corresponds to node indices (1-100) and the y-axis shows the mean FA.
    Error bars represent the standard error across subjects.
    """

    left_list = []
    right_list = []
    nodes = None

    # Loop over each subject and load the per-node FA stats if available
    for subj in subject_dirs:
        paths = get_subject_paths_wo_ses(subj)
        left_stats_file = os.path.join(paths["bun_stats_dir"], f"{track_base}_left_n100_fa_statistics.csv")
        right_stats_file = os.path.join(paths["bun_stats_dir"], f"{track_base}_right_n100_fa_statistics.csv")

        left_df = pd.read_csv(left_stats_file)
        right_df = pd.read_csv(right_stats_file)

        # Assume nodes are identical across subjects; set once
        if nodes is None:
            nodes = left_df["Node"]

        left_list.append(left_df["Mean_FA"].values)
        right_list.append(right_df["Mean_FA"].values)


    if len(left_list) == 0 or len(right_list) == 0:
        print("No subject stats files found. Skipping plot.")
        return

    left_array = np.array(left_list)  # shape: (n_subjects, 100)
    right_array = np.array(right_list)

    # Compute group-level mean and standard error (SE) across subjects
    group_mean_left = np.mean(left_array, axis=0)
    group_se_left = np.std(left_array, axis=0, ddof=1) / np.sqrt(left_array.shape[0])
    group_mean_right = np.mean(right_array, axis=0)
    group_se_right = np.std(right_array, axis=0, ddof=1) / np.sqrt(right_array.shape[0])

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.errorbar(nodes, group_mean_left, yerr=group_se_left, label=f'Left {track_base}',
                 marker='o', linestyle='-', linewidth=2, capsize=3)
    plt.errorbar(nodes, group_mean_right, yerr=group_se_right, label=f'Right {track_base}',
                 marker='s', linestyle='-', linewidth=2, capsize=3)
    plt.xlabel("Node index (1-100)")
    plt.ylabel("Mean FA Value")
    plt.title(f"Mean FA Along Nodes ({track_base} Left vs. Right) in 5 pre-TMS patients")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_paired_mean_fa_non_resampled(subject_dirs, track_base):
    """
    Gathers non-resampled mean FA values (left vs. right) for each subject
    and creates a paired (connected) dot plot. Also performs a Wilcoxon
    signed-rank test on the left vs. right arrays.

    Assumes each subject directory has a subfolder named 'bundle_stats'
    containing:
        {track_base}_left_mean_fa.csv
        {track_base}_right_mean_fa.csv
    each with at least one row and a 'Mean_FA' column.

    Example file path:
      /.../Bundles/11/bundle_stats/AF_left_mean_fa.csv

    Parameters
    ----------
    subject_dirs : list of str
        Paths to subject directories (e.g., ['/.../Bundles/11', '/.../Bundles/13', ...])
    track_base : str
        The base name of the tract (e.g., 'AF' or 'CST').

    Returns
    -------
    None
    """

    left_vals = []
    right_vals = []

    # 1) Collect left/right mean FA for each subject
    for subj in subject_dirs:
        stats_dir = os.path.join(subj, "bundle_stats")
        left_file = os.path.join(stats_dir, f"{track_base}_left_mean_fa.csv")
        right_file = os.path.join(stats_dir, f"{track_base}_right_mean_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        try:
            df_left = pd.read_csv(left_file)
            df_right = pd.read_csv(right_file)
            left_val = df_left["Mean_FA"].iloc[0]
            right_val = df_right["Mean_FA"].iloc[0]
            left_vals.append(left_val)
            right_vals.append(right_val)
        except Exception as e:
            print(f"Error reading files for subject {subj}: {e}. Skipping.")

    # 2) Check if we have any data
    if len(left_vals) == 0:
        print("No valid subjects found. Nothing to plot.")
        return

    # 3) Create the paired dot plot
    plt.figure(figsize=(6, 6))
    x_positions = [1, 2]  # Left at x=1, Right at x=2

    # Plot each subject's line from left to right
    for i in range(len(left_vals)):
        plt.plot(x_positions, [left_vals[i], right_vals[i]],
                 marker='o', color='gray', alpha=0.6, linewidth=1)

    # Plot group means as large markers
    mean_left = np.mean(left_vals)
    mean_right = np.mean(right_vals)
    plt.plot(1, mean_left, marker='o', color='blue', markersize=12,
             label=f"Mean Left = {mean_left:.3f}")
    plt.plot(2, mean_right, marker='o', color='orange', markersize=12,
             label=f"Mean Right = {mean_right:.3f}")

    plt.xlim(0.5, 2.5)
    plt.xticks([1, 2], ["Left", "Right"], fontsize=12)
    plt.ylabel("Mean FA", fontsize=12)
    plt.title(f"Non-Resampled Mean FA (Left vs. Right) - {track_base}", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_histogram_mean_fa_non_resampled(subject_dirs, track_base):
    """
    Gathers non-resampled mean FA values (left vs. right) for each subject
    and plots a histogram showing both distributions overlapped.
    """

    left_vals = []
    right_vals = []

    # Collect left/right mean FA for each subject
    for subj in subject_dirs:
        stats_dir = os.path.join(subj, "bundle_stats")
        left_file = os.path.join(stats_dir, f"{track_base}_left_mean_fa.csv")
        right_file = os.path.join(stats_dir, f"{track_base}_right_mean_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        df_left = pd.read_csv(left_file)
        df_right = pd.read_csv(right_file)

        left_val = df_left["Mean_FA"].iloc[0]
        right_val = df_right["Mean_FA"].iloc[0]

        left_vals.append(left_val)
        right_vals.append(right_val)

    # Check if we have any data
    if len(left_vals) == 0:
        print("No valid subjects found. Nothing to plot.")
        return
    print("Number of valid subjects:", len(left_vals))
    print("Left FA values:", left_vals)
    print("Right FA values:", right_vals)

    # Plot overlapping histograms
    plt.figure(figsize=(7, 5))
    plt.hist(left_vals, bins=10, alpha=0.5, color='blue', label='Left')
    plt.hist(right_vals, bins=10, alpha=0.5, color='orange', label='Right')
    plt.xlabel("Mean FA")
    plt.ylabel("Frequency")
    plt.title(f"Histogram of Non-Resampled Mean FA for {track_base} (Left vs. Right)")
    plt.legend()

    plt.show()
