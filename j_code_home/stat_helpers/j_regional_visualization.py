import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from helpers.j_helpers import get_subject_paths


def plot_mean_metrics_along_nodes(subject_dirs, track_base, metrics, alpha=1, save_figure=False):
    """
    Aggregates node-wise mean metric values across subjects for the left and right
    resampled tracts and creates multiple subplots:
      - One subplot per non-peaks metric (e.g., FA, ADC).
      - Two subplots for peaks (original vs. absolute).

    Example layout if metrics = [FA, ADC, peaks] (2x2):
      ┌────────────┬────────────┐
      │  FA Plot   │  ADC Plot  │
      ├────────────┼────────────┤
      │  Orig PEAKS│  Abs PEAKS │
      └────────────┴────────────┘

    Each subplot shows node index (1-100) on the x-axis and the mean metric on the y-axis,
    with error bars for the standard error across subjects.
    """
    # Separate 'peaks' from the other metrics
    normal_metrics = [m for m in metrics if m.lower() != "peaks"]
    has_peaks = any(m.lower() == "peaks" for m in metrics)

    # Each normal metric gets 1 subplot; peaks gets 2 subplots
    n_subplots = len(normal_metrics) + (2 if has_peaks else 0)
    if n_subplots == 0:
        print("No metrics provided. Exiting.")
        return

    # Determine a suitable grid layout (e.g., 2x2 if we have 4 subplots)
    nrows = math.ceil(math.sqrt(n_subplots))
    ncols = math.ceil(n_subplots / nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    # Flatten or wrap the axes array for easy iteration
    if n_subplots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    i = 0  # index to track which subplot we're filling
    sessions = ['ses_pre']

    ############################################################################
    # 1) Plot each normal metric (e.g., FA, ADC) in its own subplot
    ############################################################################
    for metric in normal_metrics:
        ax = axes[i]
        metric_lower = metric.lower()
        ext = "csv"  # we assume stats are always CSV for FA/ADC
        left_list = []
        right_list = []
        nodes = None

        # Gather node-wise stats across subjects
        for subj_dir in subject_dirs:
            for session in sessions:
                session_dir = os.path.join(subj_dir, session)
                paths = get_subject_paths(session_dir)
                left_stats_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_left_n100_{metric_lower}_statistics.{ext}"
                )
                right_stats_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_right_n100_{metric_lower}_statistics.{ext}"
                )
                if not (os.path.exists(left_stats_file) and os.path.exists(right_stats_file)):
                    print(f"Missing stats file(s) for {metric} in {session_dir}. Skipping.")
                    continue

                left_df = pd.read_csv(left_stats_file)
                right_df = pd.read_csv(right_stats_file)

                if nodes is None:
                    nodes = left_df["Node"]  # assume same nodes across subjects
                left_list.append(left_df[f"Mean_{metric.upper()}"].values)
                right_list.append(right_df[f"Mean_{metric.upper()}"].values)

        if len(left_list) == 0 or len(right_list) == 0:
            print(f"No valid stats files found for {metric}. Skipping subplot.")
            i += 1
            continue

        # Convert to numpy arrays and compute group means + standard error
        left_array = np.array(left_list)
        right_array = np.array(right_list)
        group_mean_left = np.mean(left_array, axis=0)
        group_se_left = np.std(left_array, axis=0, ddof=1) / np.sqrt(left_array.shape[0])
        group_mean_right = np.mean(right_array, axis=0)
        group_se_right = np.std(right_array, axis=0, ddof=1) / np.sqrt(right_array.shape[0])

        # Plot
        ax.errorbar(nodes, group_mean_left, yerr=group_se_left,
                    label=f'Left {track_base}', marker='o', linestyle='-', linewidth=2, capsize=3, alpha=alpha)
        ax.errorbar(nodes, group_mean_right, yerr=group_se_right,
                    label=f'Right {track_base}', marker='s', linestyle='-', linewidth=2, capsize=3, alpha=alpha)
        ax.set_xlabel("Node index (1-100)")
        ax.set_ylabel(f"Mean {metric.upper()} Value")
        ax.set_title(f"Mean {metric.upper()} Along Nodes ({track_base} Left vs. Right)")
        ax.legend()
        ax.grid(True)
        i += 1

    ############################################################################
    # 2) If we have peaks, create two subplots: one for original, one for abs
    ############################################################################
    if has_peaks:
        # We'll collect original and absolute data
        left_orig_list, right_orig_list = [], []
        left_abs_list, right_abs_list = [], []
        nodes = None  # same nodes for both

        for subj_dir in subject_dirs:
            for session in sessions:
                session_dir = os.path.join(subj_dir, session)
                paths = get_subject_paths(session_dir)
                # Original peaks stats
                left_orig_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_left_n100_peaks_statistics.csv"
                )
                right_orig_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_right_n100_peaks_statistics.csv"
                )
                # Absolute peaks stats
                left_abs_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_left_n100_peaks_abs_statistics.csv"
                )
                right_abs_file = os.path.join(
                    paths["bun_stats_dir"], f"{track_base}_right_n100_peaks_abs_statistics.csv"
                )

                # Make sure all four files exist
                if not all(os.path.exists(f) for f in [left_orig_file, right_orig_file,
                                                       left_abs_file, right_abs_file]):
                    print(f"Missing peaks stats file(s) for subject {session_dir}. Skipping.")
                    continue

                left_orig_df = pd.read_csv(left_orig_file)
                right_orig_df = pd.read_csv(right_orig_file)
                left_abs_df = pd.read_csv(left_abs_file)
                right_abs_df = pd.read_csv(right_abs_file)

                if nodes is None:
                    nodes = left_orig_df["Node"]  # assume same nodes across subjects
                left_orig_list.append(left_orig_df["Mean_PEAKS"].values)
                right_orig_list.append(right_orig_df["Mean_PEAKS"].values)
                left_abs_list.append(left_abs_df["Mean_PEAKS"].values)
                right_abs_list.append(right_abs_df["Mean_PEAKS"].values)

        if (not left_orig_list or not right_orig_list or
                not left_abs_list or not right_abs_list):
            print("No valid peaks stats files found. Skipping peaks subplots.")
            return

        # Convert to arrays
        left_orig_array = np.array(left_orig_list)
        right_orig_array = np.array(right_orig_list)
        left_abs_array = np.array(left_abs_list)
        right_abs_array = np.array(right_abs_list)

        # Compute group-level means and SE for original
        group_mean_left_orig = np.mean(left_orig_array, axis=0)
        group_se_left_orig = np.std(left_orig_array, axis=0, ddof=1) / np.sqrt(left_orig_array.shape[0])
        group_mean_right_orig = np.mean(right_orig_array, axis=0)
        group_se_right_orig = np.std(right_orig_array, axis=0, ddof=1) / np.sqrt(right_orig_array.shape[0])

        # Compute group-level means and SE for absolute
        group_mean_left_abs = np.mean(left_abs_array, axis=0)
        group_se_left_abs = np.std(left_abs_array, axis=0, ddof=1) / np.sqrt(left_abs_array.shape[0])
        group_mean_right_abs = np.mean(right_abs_array, axis=0)
        group_se_right_abs = np.std(right_abs_array, axis=0, ddof=1) / np.sqrt(right_abs_array.shape[0])

        # (a) Subplot for ORIGINAL peaks
        ax_orig = axes[i]
        ax_orig.errorbar(nodes, group_mean_left_orig, yerr=group_se_left_orig,
                         label=f'Left {track_base}', marker='o', linestyle='-', linewidth=2, capsize=3, alpha=alpha)
        ax_orig.errorbar(nodes, group_mean_right_orig, yerr=group_se_right_orig,
                         label=f'Right {track_base}', marker='s', linestyle='-', linewidth=2, capsize=3, alpha=alpha)
        ax_orig.set_xlabel("Node index (1-100)")
        ax_orig.set_ylabel("Mean peaks value")
        ax_orig.set_title(f"Mean peaks ({track_base} Left vs. Right)")
        ax_orig.legend()
        ax_orig.grid(True)
        i += 1

        # (b) Subplot for ABS peaks
        ax_abs = axes[i]
        ax_abs.errorbar(nodes, group_mean_left_abs, yerr=group_se_left_abs,
                        label=f'Left {track_base} abs', marker='^', linestyle='--', linewidth=2, capsize=3, alpha=alpha)
        ax_abs.errorbar(nodes, group_mean_right_abs, yerr=group_se_right_abs,
                        label=f'Right {track_base} abs', marker='v', linestyle='--', linewidth=2, capsize=3, alpha=alpha)
        ax_abs.set_xlabel("Node index (1-100)")
        ax_abs.set_ylabel("Mean peaks value")
        ax_abs.set_title(f"Mean absolute peaks ({track_base} Left vs. Right)")
        ax_abs.legend()
        ax_abs.grid(True)
        i += 1

    # Hide any unused subplots (if any)
    for j in range(i, len(axes)):
        axes[j].axis("off")

    title_text = "Interhemispheric Comparison of Resampled AF Metrics (FA, ADC, and Peaks) in 5 Healthy Subjects"
    plt.suptitle(title_text, fontsize=14)
    plt.tight_layout()
    plt.show()

    if save_figure:
        save_dir = "/Users/nikitakaruzin/Desktop/Research/Picht/Lab Report/figs/bun_figs"
        os.makedirs(save_dir, exist_ok=True)  # ensure directory exists

        # Create a safe filename from the title by replacing spaces with underscores
        safe_filename = title_text.replace(' ', '_') + ".png"
        save_path = os.path.join(save_dir, safe_filename)

        # Save with higher resolution (dpi=300) and tight bounding box
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")


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
