import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from helpers.j_helpers import (run_cmd, create_tractseg_file,
                               get_subject_dirs, get_subject_paths, get_args, logs)

###############################################
# 1. Process Mean Metric for non-resampled tracts #
###############################################
def process_track_mean(paths, subject_dir, track_base, metric):
    """
    Opens the metric file for the left and right tracts (e.g., AF or CST),
    computes the mean value for each tract, and saves the output as two CSV files.

    Expected file paths (within subject_dir):
      - ses_pre/along_tract/<metric>/<track_base>_left_<metric_lower>.<ext>
      - ses_pre/along_tract/<metric>/<track_base>_right_<metric_lower>.<ext>

    For "peaks", the file extension is .txt; for other metrics, it is .csv.

    Output files (saved in paths["bun_stats_dir"]):
      - <track_base>_left_mean_<metric_lower>.csv
      - <track_base>_right_mean_<metric_lower>.csv
    """
    metric_lower = metric.lower()
    ext = "txt" if metric_lower == "peaks" else "csv"
    metric_dir = os.path.join(subject_dir, "along_tract", metric)
    if not os.path.exists(metric_dir):
        print(f"{metric} directory not found at {metric_dir}.")
        return

    os.makedirs(paths["bun_stats_dir"], exist_ok=True)

    # Loop over left and right tracts
    for side in ["left", "right"]:
        tract = f"{track_base}_{side}"
        file_path = os.path.join(metric_dir, f"{tract}_{metric_lower}.{ext}")
        if not os.path.exists(file_path):
            print(f"File {file_path} not found. Skipping {tract}.")
            continue

        all_values = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.replace('\x00', '').strip()
                if not line or line.startswith('#'):
                    continue
                # Convert each whitespace‐separated number to float
                values = [float(x) for x in line.split()]
                all_values.extend(values)

        mean_val = np.mean(all_values)
        df = pd.DataFrame({'Tract': [tract], f"Mean_{metric.upper()}": [mean_val]})
        output_filename = f"{tract}_mean_{metric_lower}.csv"
        output_file = os.path.join(paths["bun_stats_dir"], output_filename)
        df.to_csv(output_file, index=False)
        print(f"Mean {metric} for {tract} saved to {output_file}")

##########################################
# 2. Process Metric Stats for resampled data #
##########################################
def process_metric_stats(paths, input_file, metric):
    """
    Reads a file containing metric values for streamlines along nodes,
    computes node-wise summary statistics, and saves the results into a CSV file under paths["bun_stats_dir"].
    """
    data = pd.read_csv(input_file, header=None, comment='#', sep=r'\s+')
    data = data.astype(float)

    mean_per_node = data.mean(axis=0)
    std_per_node = data.std(axis=0)

    base = os.path.basename(input_file)
    name, _ = os.path.splitext(base)
    output_filename = f"{name}_statistics.csv"
    os.makedirs(paths["bun_stats_dir"], exist_ok=True)
    output_file = os.path.join(paths["bun_stats_dir"], output_filename)

    stats_df = pd.DataFrame({
        "Node": mean_per_node.index + 1,
        f"Mean_{metric.upper()}": mean_per_node.values,
        f"Std_{metric.upper()}": std_per_node.values
    })
    stats_df.to_csv(output_file, index=False)
    print(f"Statistics saved to {output_file}")
    return {"mean_per_node": mean_per_node, "std_per_node": std_per_node, "data": data}

##########################################
# 3. Run Metric Processing for resampled data #
##########################################
def run_metric_processing(paths, subj_dir, track_base, metric):
    """
    Processes metric files for specified resampled tracts for the current subject/session.

    Expected metric file names (in along_tract/<metric> folder):
      - <track_base>_left_n100_<metric_lower>.<ext>
      - <track_base>_right_n100_<metric_lower>.<ext>
    """
    metric_lower = metric.lower()
    ext = "txt" if metric_lower == "peaks" else "csv"
    session_dir = subj_dir
    metric_dir = os.path.join(session_dir, "along_tract", metric)
    if not os.path.exists(metric_dir):
        print(f"{metric} directory {metric_dir} not found. Skipping processing for current session.")
        return

    for side in ["left", "right"]:
        tract = f"{track_base}_{side}_n100"
        file_path = os.path.join(metric_dir, f"{tract}_{metric_lower}.{ext}")
        if os.path.exists(file_path):
            print(f"Processing {file_path}...")
            process_metric_stats(paths, file_path, metric)
        else:
            print(f"File {file_path} not found in {metric_dir}. Skipping.")

#################################################
# 4. Group-wise paired t-test for resampled data #
#################################################
def groupwise_track_ttest(root, subject_dirs, track_base, metric):
    """
    Performs a node-wise paired t-test comparing metric values for the left and right
    tracts (e.g., AF or CST) across subjects (using only the ses_pre session).

    Expected files for each subject:
      - <subject_dir>/ses_pre/along_tract/<metric>/<track_base>_left_n100_<metric_lower>.<ext>
      - <subject_dir>/ses_pre/along_tract/<metric>/<track_base>_right_n100_<metric_lower>.<ext>

    The results are saved to:
      <root>/group_analysis/groupwise_<track_base>_left_vs_right_ttest_<metric_lower>.csv
    """
    metric_lower = metric.lower()
    ext = "txt" if metric_lower == "peaks" else "csv"
    left_list = []
    right_list = []

    for subj_dir in subject_dirs:
        sessions = ['ses_pre']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            metric_dir = os.path.join(session_dir, "along_tract", metric)
            left_file = os.path.join(metric_dir, f"{track_base}_left_n100_{metric_lower}.{ext}")
            right_file = os.path.join(metric_dir, f"{track_base}_right_n100_{metric_lower}.{ext}")

            if not os.path.exists(left_file) or not os.path.exists(right_file):
                print(f"Missing file(s) for subject {session_dir}. Skipping.")
                continue

            try:
                left_df = pd.read_csv(left_file, comment='#', header=None, sep='\s+')
                right_df = pd.read_csv(right_file, comment='#', header=None, sep='\s+')
            except Exception as e:
                print(f"Error reading files for subject {session_dir}: {e}. Skipping.")
                continue

            if left_df.shape[1] != 100 or right_df.shape[1] != 100:
                print(f"Subject {session_dir} does not have 100 nodes. Skipping. "
                      f"Left shape: {left_df.shape}, Right shape: {right_df.shape}")
                continue

            left_mean = left_df.mean(axis=0).values
            right_mean = right_df.mean(axis=0).values
            left_list.append(left_mean)
            right_list.append(right_mean)

    if len(left_list) == 0:
        print("No valid subjects found for group-wise t-test.")
        return

    left_array = np.array(left_list)
    right_array = np.array(right_list)
    n_subjects, n_nodes = left_array.shape
    print(f"Performing paired t-tests across {n_subjects} subjects and {n_nodes} nodes for {track_base} with {metric}.")

    t_stats = []
    p_values = []
    for node in range(n_nodes):
        t, p = stats.ttest_rel(left_array[:, node], right_array[:, node])
        t_stats.append(t)
        p_values.append(p)

    _, p_values_fdr, _, _ = multipletests(p_values, method='fdr_bh')
    results_df = pd.DataFrame({
        "Node": np.arange(1, n_nodes + 1),
        "t_stat": t_stats,
        "p_value": p_values,
        "p_value_fdr": p_values_fdr
    })

    output_dir = os.path.join(root, "group_analysis")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"groupwise_{track_base}_left_vs_right_ttest_{metric_lower}.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Group-wise paired t-test results saved to {output_file}")

###########################################################
# 5. Group-wise Wilcoxon Signed-Ranks test for resampled data #
###########################################################
def groupwise_track_wilcoxon(root, subject_dirs, track_base, metric):
    """
    Performs a node-wise Wilcoxon Signed-Ranks Test comparing metric values for the left and right
    tracts (e.g., AF or CST) across subjects.

    Expected files for each subject:
      - <subject_dir>/ses_pre/along_tract/<metric>/<track_base>_left_n100_<metric_lower>.<ext>
      - <subject_dir>/ses_pre/along_tract/<metric>/<track_base>_right_n100_<metric_lower>.<ext>

    The results are saved to:
      <root>/group_analysis/groupwise_<track_base>_left_vs_right_wilcoxon_<metric_lower>.csv
    """
    metric_lower = metric.lower()
    ext = "txt" if metric_lower == "peaks" else "csv"
    left_list = []
    right_list = []

    for subj_dir in subject_dirs:
        sessions = ['ses_pre']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            metric_dir = os.path.join(session_dir, "along_tract", metric)
            left_file = os.path.join(metric_dir, f"{track_base}_left_n100_{metric_lower}.{ext}")
            right_file = os.path.join(metric_dir, f"{track_base}_right_n100_{metric_lower}.{ext}")

            if not os.path.exists(left_file) or not os.path.exists(right_file):
                print(f"Missing file(s) for subject {session_dir}. Skipping.")
                continue

            left_df = pd.read_csv(left_file, comment='#', header=None, sep='\s+')
            right_df = pd.read_csv(right_file, comment='#', header=None, sep='\s+')

            if left_df.shape[1] != 100 or right_df.shape[1] != 100:
                print(f"Subject {session_dir} does not have 100 nodes. Skipping. "
                      f"Left shape: {left_df.shape}, Right shape: {right_df.shape}")
                continue

            left_mean = left_df.mean(axis=0).values
            right_mean = right_df.mean(axis=0).values
            left_list.append(left_mean)
            right_list.append(right_mean)

    if len(left_list) == 0:
        print("No valid subjects found for group-wise Wilcoxon test.")
        return

    left_array = np.array(left_list)
    right_array = np.array(right_list)
    n_subjects, n_nodes = left_array.shape
    print(f"Performing Wilcoxon signed-ranks tests across {n_subjects} subjects and {n_nodes} nodes (RESAMPLED)"
          f" for {track_base} with {metric}.")

    stat_vals = []
    p_values = []
    for node in range(n_nodes):
        stat, p = stats.wilcoxon(left_array[:, node], right_array[:, node])
        stat_vals.append(stat)
        p_values.append(p)

    _, p_values_fdr, _, _ = multipletests(p_values, method='fdr_bh')
    results_df = pd.DataFrame({
        "Node": np.arange(1, n_nodes + 1),
        "statistic": stat_vals,
        "p_value": p_values,
        "p_value_fdr": p_values_fdr
    })

    output_dir = os.path.join(root, "group_analysis", "resampled")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"groupwise_{track_base}_left_vs_right_wilcoxon_{metric_lower}_resampled.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Group-wise Wilcoxon test results saved to {output_file}")

    alpha = 0.05
    df = pd.read_csv(output_file)
    n_uncorrected = (df["p_value"] < alpha).sum()
    n_fdr_corrected = (df["p_value_fdr"] < alpha).sum()
    min_p_value = df["p_value"].min()
    max_p_value = df["p_value"].max()

    print(f"Number of nodes with uncorrected p < {alpha}: {n_uncorrected}")
    print(f"Number of nodes with FDR-corrected p < {alpha}: {n_fdr_corrected}")
    print(f"Minimum uncorrected p-value: {min_p_value}")
    print(f"Maximum uncorrected p-value: {max_p_value}")

###############################################
# 6. Group-wise Wilcoxon test for non-resampled data #
###############################################
def groupwise_track_wilcoxon_non_resampled(root, subject_dirs, track_base, metric):
    """
    Performs a Wilcoxon Signed-Ranks Test on the non-resampled (single) mean metric
    values of the left vs. right tract (e.g., AF or CST) across subjects.

    For each subject directory, it looks for the following files in:
        <subject_dir>/ses_pre/bundle_stats/<track_base>_left_mean_<metric_lower>.csv
        <subject_dir>/ses_pre/bundle_stats/<track_base>_right_mean_<metric_lower>.csv

    Each CSV should contain at least one row with a column "Mean_<METRIC>".
    The function:
      1. Reads each subject's left/right mean metric,
      2. Performs a paired Wilcoxon test (left vs. right),
      3. Saves a one-row CSV file with the statistic, p-value, number of subjects,
         and group means to:
         <root>/group_analysis/groupwise_<track_base>_left_vs_right_wilcoxon_non_resampled_<metric_lower>.csv
    """
    metric_lower = metric.lower()
    left_vals = []
    right_vals = []

    for subj_dir in subject_dirs:
        sessions = ['ses_pre']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            stats_dir = os.path.join(session_dir, "bundle_stats")
            left_file = os.path.join(stats_dir, f"{track_base}_left_mean_{metric_lower}.csv")
            right_file = os.path.join(stats_dir, f"{track_base}_right_mean_{metric_lower}.csv")

            if not os.path.exists(left_file) or not os.path.exists(right_file):
                print(f"Missing file(s) for subject {session_dir}. Skipping.")
                continue

            try:
                df_left = pd.read_csv(left_file)
                df_right = pd.read_csv(right_file)
                # Use the metric-specific column name if present; otherwise fallback
                left_val = df_left[f"Mean_{metric.upper()}"].iloc[0] if f"Mean_{metric.upper()}" in df_left.columns else df_left["Mean_FA"].iloc[0]
                right_val = df_right[f"Mean_{metric.upper()}"].iloc[0] if f"Mean_{metric.upper()}" in df_right.columns else df_right["Mean_FA"].iloc[0]
                left_vals.append(left_val)
                right_vals.append(right_val)
            except Exception as e:
                print(f"Error reading files for subject {session_dir}: {e}. Skipping.")
                continue

    if len(left_vals) == 0:
        print("No valid subjects found for group-wise Wilcoxon test (non-resampled).")
        return

    stat, p_value = stats.wilcoxon(left_vals, right_vals)
    group_mean_left = np.mean(left_vals)
    group_mean_right = np.mean(right_vals)

    output_dir = os.path.join(root, "group_analysis", "non_resampled")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(
        output_dir, f"groupwise_{track_base}_left_vs_right_wilcoxon_{metric_lower}_non_resampled.csv"
    )

    results_df = pd.DataFrame({
        "Wilcoxon_Statistic": [stat],
        "p_value": [p_value],
        "N_subjects": [len(left_vals)],
        "Group_mean_left": [group_mean_left],
        "Group_mean_right": [group_mean_right]
    })

    results_df.to_csv(output_file, index=False)
    print(f"Group-wise non-resampled Wilcoxon results saved to {output_file}")
    print(f"Wilcoxon statistic: {stat}, p-value: {p_value}")
    print(f"Mean left: {group_mean_left}, Mean right: {group_mean_right}")

###############################################
# 7. Plot Tractometry Results (unchanged) #
###############################################
def plot_tractometry_results(root, bundles="AF_left AF_right"):
    """
    This function first creates the necessary TractSeg file with tractometry path, bundles to use,
    and a 3D plot option. It then creates a png with the txt options.
    """
    tractseg_tractometry_dir = os.path.join(root, "group_analysis", "tractseg_tractometry")
    subject_txt_output = os.path.join(tractseg_tractometry_dir, "subjects.txt")
    subjects_txt = create_tractseg_file(
        root,
        tractometry_path=f"{tractseg_tractometry_dir}/SUBJECT_ID/Tractometry.csv",
        bundles=bundles,
        plot3d=tractseg_tractometry_dir,
        output_file=subject_txt_output)
    tractometry_png = os.path.join(root, "group_analysis", "tractseg_tractometry", "tractometry_result.png")
    run_cmd([
        "plot_tractometry_results",
        "-i", subjects_txt,
        "-o", tractometry_png,
        "--mc"
    ])

###############################################
# 8. Get Node-wise Metric Range Across Subjects #
###############################################
def get_nodewise_metric_range(subject_dirs, track_base, metric):
    """
    Gathers all node-wise metric values for the left and right tracts (e.g.,
    <track_base>_left_n100_<metric_lower>_statistics.csv) across all subjects, then computes the global min and max.

    Expected file naming convention for each subject (under ses_pre/bundle_stats):
      - <track_base>_left_n100_<metric_lower>_statistics.csv
      - <track_base>_right_n100_<metric_lower>_statistics.csv

    Each CSV should have columns: ["Node", "Mean_<METRIC>", "Std_<METRIC>"].
    """
    metric_lower = metric.lower()
    all_left_values = []
    all_right_values = []

    for subj_dir in subject_dirs:
        sessions = ['ses_pre']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            stats_dir = os.path.join(session_dir, "bundle_stats")
            left_stats_file = os.path.join(stats_dir, f"{track_base}_left_n100_{metric_lower}_statistics.csv")
            right_stats_file = os.path.join(stats_dir, f"{track_base}_right_n100_{metric_lower}_statistics.csv")

            if not os.path.exists(left_stats_file) or not os.path.exists(right_stats_file):
                print(f"Missing stats file(s) in {stats_dir}. Skipping subject {session_dir}.")
                continue

            try:
                df_left = pd.read_csv(left_stats_file)
                df_right = pd.read_csv(right_stats_file)
                all_left_values.extend(df_left[f"Mean_{metric.upper()}"].values)
                all_right_values.extend(df_right[f"Mean_{metric.upper()}"].values)
            except Exception as e:
                print(f"Error reading stats for subject {session_dir}: {e}. Skipping.")

    if len(all_left_values) == 0 or len(all_right_values) == 0:
        print("No valid node-wise data found across subjects.")
        return None

    left_min, left_max = min(all_left_values), max(all_left_values)
    right_min, right_max = min(all_right_values), max(all_right_values)

    print(f"\n--- {track_base} Node-wise {metric.upper()} Range Across All Subjects ---")
    print(f"Left {track_base}:  min = {left_min:.4f}, max = {left_max:.4f}")
    print(f"Right {track_base}: min = {right_min:.4f}, max = {right_max:.4f}\n")

    return left_min, left_max, right_min, right_max

###############################################
# 9. Report Non-resampled Metric Stats #
###############################################
def report_nonresampled_metric_stats(subject_dirs, track_base, metric):
    """
    Gathers non-resampled mean metric values for left vs. right tracts across subjects,
    computes group mean ± SD, runs a Wilcoxon signed-rank test, and prints a summary.

    Expected files for each subject (under ses_pre/bundle_stats):
      - <track_base>_left_mean_<metric_lower>.csv
      - <track_base>_right_mean_<metric_lower>.csv

    Each CSV must have a column "Mean_<METRIC>" in the first row.
    """
    metric_lower = metric.lower()
    left_vals = []
    right_vals = []

    sessions = ['ses_pre']
    for subj_dir in subject_dirs:
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            stats_dir = os.path.join(session_dir, "bundle_stats")
            left_file = os.path.join(stats_dir, f"{track_base}_left_mean_{metric_lower}.csv")
            right_file = os.path.join(stats_dir, f"{track_base}_right_mean_{metric_lower}.csv")

            if not os.path.exists(left_file) or not os.path.exists(right_file):
                print(f"Missing file(s) for subject {session_dir}. Skipping.")
                continue

            try:
                df_left = pd.read_csv(left_file)
                df_right = pd.read_csv(right_file)
                left_val = df_left[f"Mean_{metric.upper()}"].iloc[0] if f"Mean_{metric.upper()}" in df_left.columns else df_left["Mean_FA"].iloc[0]
                right_val = df_right[f"Mean_{metric.upper()}"].iloc[0] if f"Mean_{metric.upper()}" in df_right.columns else df_right["Mean_FA"].iloc[0]
                left_vals.append(left_val)
                right_vals.append(right_val)
            except Exception as e:
                print(f"Error reading files for subject {session_dir}: {e}. Skipping.")

    if len(left_vals) == 0:
        print("No valid subjects found for non-resampled metric. Cannot compute summary.")
        return

    mean_left = np.mean(left_vals)
    std_left = np.std(left_vals, ddof=1)
    mean_right = np.mean(right_vals)
    std_right = np.std(right_vals, ddof=1)

    stat, p_value = stats.wilcoxon(left_vals, right_vals)

    print(
        f"When averaging along the non-resampled tract, the mean left {track_base} {metric.upper()} was "
        f"{mean_left:.3f} (± {std_left:.3f}) compared with {mean_right:.3f} (± {std_right:.3f}) "
        f"on the right (Wilcoxon signed-rank p={p_value:.3f})."
    )


def convert_peaks_to_absolute(paths, session_dir, track_base):
    """
    Reads the n100 resampled peaks files (e.g., <track_base>_left_n100_peaks.txt),
    converts the values to absolute values, and saves them as new *_abs.txt files.
    Then it calls process_metric_stats on these *_abs.txt files to produce
    *_abs_statistics.csv in paths["bun_stats_dir"].

    Expected input files (within session_dir):
      - along_tract/peaks/<track_base>_left_n100_peaks.txt
      - along_tract/peaks/<track_base>_right_n100_peaks.txt

    Output files:
      - <track_base>_left_n100_peaks_abs.txt
      - <track_base>_right_n100_peaks_abs.txt
      - plus their corresponding *_abs_statistics.csv files (via process_metric_stats).
    """
    import os
    from stat_helpers.j_regional_analysis import process_metric_stats  # or wherever process_metric_stats is defined

    metric = "peaks"
    metric_dir = os.path.join(session_dir, "along_tract", metric)
    if not os.path.exists(metric_dir):
        print(f"Peaks directory not found at {metric_dir}.")
        return

    # Only handle the n100 files (left/right)
    for side in ["left", "right"]:
        input_filename = f"{track_base}_{side}_n100_{metric}.txt"
        input_file = os.path.join(metric_dir, input_filename)
        if not os.path.exists(input_file):
            print(f"File {input_file} not found. Skipping.")
            continue

        with open(input_file, 'r') as f:
            lines = f.readlines()

        abs_lines = []
        for line in lines:
            line_stripped = line.strip()
            # Keep empty or commented lines as-is
            if not line_stripped or line_stripped.startswith("#"):
                abs_lines.append(line)
                continue
            try:
                # Convert each whitespace-separated number to abs value
                values = [abs(float(x)) for x in line_stripped.split()]
                abs_line = " ".join(str(val) for val in values) + "\n"
                abs_lines.append(abs_line)
            except Exception as e:
                print(f"Error processing line '{line_stripped}': {e}. Keeping original line.")
                abs_lines.append(line)

        # Write out the new _abs.txt file
        output_filename = f"{track_base}_{side}_n100_{metric}_abs.txt"
        output_file = os.path.join(metric_dir, output_filename)
        with open(output_file, 'w') as f:
            f.writelines(abs_lines)

        print(f"Absolute peak values for {output_filename} created.")

        # Immediately run process_metric_stats on the new file to get *_abs_statistics.csv
        process_metric_stats(paths, input_file=output_file, metric="peaks")
