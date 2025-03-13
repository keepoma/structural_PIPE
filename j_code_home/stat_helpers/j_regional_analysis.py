import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from helpers.j_helpers import (run_cmd, create_tractseg_file,
                               get_subject_dirs, get_subject_paths_wo_ses, get_args, logs)

###############################################
# 1. Process Mean FA for non-resampled tracts #
###############################################
def process_track_mean_fa(paths, subject_dir, track_base):
    """
    Opens the FA CSV files for the left and right tracts (e.g., AF or CST),
    computes the mean FA for each tract, and saves the output as two CSV files.

    Expected file paths (within subject_dir):
      - ses_pre/along_tract/FA/<track_base>_left_fa.csv
      - ses_pre/along_tract/FA/<track_base>_right_fa.csv

    Output files (saved in paths["bun_stats_dir"]):
      - <track_base>_left_mean_fa.csv
      - <track_base>_right_mean_fa.csv
    """
    fa_dir = os.path.join(subject_dir, "along_tract", "FA")
    if not os.path.exists(fa_dir):
        print(f"FA directory not found at {fa_dir}.")
        return

    os.makedirs(paths["bun_stats_dir"], exist_ok=True)

    tracts = [f"{track_base}_left", f"{track_base}_right"]

    for tract in tracts:
        file_path = os.path.join(fa_dir, f"{tract}_fa.csv")
        if not os.path.exists(file_path):
            print(f"File {file_path} not found. Skipping {tract}.")
            continue

        all_values = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.replace('\x00', '').strip()
                if not line or line.startswith('#'):
                    continue
                values = [float(x) for x in line.split()]
                all_values.extend(values)

        mean_fa = np.mean(all_values)
        df = pd.DataFrame({'Tract': [tract], 'Mean_FA': [mean_fa]})
        output_filename = f"{tract}_mean_fa.csv"
        output_file = os.path.join(paths["bun_stats_dir"], output_filename)
        df.to_csv(output_file, index=False)
        print(f"Mean FA for {tract} saved to {output_file}")

##########################################
# 2. Process FA Stats for resampled data #
##########################################
def process_fa_stats(paths, input_file):
    """
    Reads a file containing FA values for streamlines along nodes,
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
        "Mean_FA": mean_per_node.values,
        "Std_FA": std_per_node.values
    })
    stats_df.to_csv(output_file, index=False)
    print(f"Statistics saved to {output_file}")
    return {"mean_per_node": mean_per_node, "std_per_node": std_per_node, "data": data}

##########################################
# 3. Run FA processing for resampled data #
##########################################
def run_fa_processing(paths, subj_dir, track_base):
    """
    Processes FA files for specified resampled tracts for the current subject/session.

    Expected FA file names (in along_tract/FA folder):
      - <track_base>_left_n100_fa.csv
      - <track_base>_right_n100_fa.csv
    """
    session_dir = subj_dir
    fa_dir = os.path.join(session_dir, "along_tract", "FA")
    if not os.path.exists(fa_dir):
        print(f"FA directory {fa_dir} not found. Skipping processing for current session.")
        return

    tract_names = [f"{track_base}_left_n100", f"{track_base}_right_n100"]
    for tract in tract_names:
        fa_csv = os.path.join(fa_dir, f"{tract}_fa.csv")
        if os.path.exists(fa_csv):
            print(f"Processing {fa_csv}...")
            process_fa_stats(paths, fa_csv)
        else:
            print(f"File {fa_csv} not found in {fa_dir}. Skipping.")

#################################################
# 4. Group-wise paired t-test for resampled data #
#################################################
def groupwise_track_ttest(root, subject_dirs, track_base):
    """
    Performs a node-wise paired t-test comparing FA values for the left and right
    tracts (e.g., AF or CST) across subjects (using only the ses_pre session).

    Expected files for each subject:
      - <subject_dir>/ses_pre/along_tract/FA/<track_base>_left_n100_fa.csv
      - <subject_dir>/ses_pre/along_tract/FA/<track_base>_right_n100_fa.csv

    The results are saved to:
      <root>/group_analysis/groupwise_<track_base>_left_vs_right_ttest.csv
    """
    left_list = []
    right_list = []

    for subj in subject_dirs:
        fa_dir = os.path.join(subj, "along_tract", "FA")
        left_file = os.path.join(fa_dir, f"{track_base}_left_n100_fa.csv")
        right_file = os.path.join(fa_dir, f"{track_base}_right_n100_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        try:
            left_df = pd.read_csv(left_file, comment='#', header=None, sep='\s+')
            right_df = pd.read_csv(right_file, comment='#', header=None, sep='\s+')
        except Exception as e:
            print(f"Error reading files for subject {subj}: {e}. Skipping.")
            continue

        if left_df.shape[1] != 100 or right_df.shape[1] != 100:
            print(f"Subject {subj} does not have 100 nodes. Skipping. "
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
    print(f"Performing paired t-tests across {n_subjects} subjects and {n_nodes} nodes for {track_base}.")

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
    output_file = os.path.join(output_dir, f"groupwise_{track_base}_left_vs_right_ttest.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Group-wise paired t-test results saved to {output_file}")

###########################################################
# 5. Group-wise Wilcoxon Signed-Ranks test for resampled data #
###########################################################
def groupwise_track_wilcoxon(root, subject_dirs, track_base):
    """
    Performs a node-wise Wilcoxon Signed-Ranks Test comparing FA values for the left and right
    tracts (e.g., AF or CST) across subjects.

    Expected files for each subject:
      - <subject_dir>/ses_pre/along_tract/FA/<track_base>_left_n100_fa.csv
      - <subject_dir>/ses_pre/along_tract/FA/<track_base>_right_n100_fa.csv

    The results are saved to:
      <root>/group_analysis/groupwise_<track_base>_left_vs_right_wilcoxon.csv
    """
    left_list = []
    right_list = []

    for subj in subject_dirs:
        fa_dir = os.path.join(subj, "along_tract", "FA")
        left_file = os.path.join(fa_dir, f"{track_base}_left_n100_fa.csv")
        right_file = os.path.join(fa_dir, f"{track_base}_right_n100_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        left_df = pd.read_csv(left_file, comment='#', header=None, sep='\s+')
        right_df = pd.read_csv(right_file, comment='#', header=None, sep='\s+')

        if left_df.shape[1] != 100 or right_df.shape[1] != 100:
            print(f"Subject {subj} does not have 100 nodes. Skipping. "
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
    print(f"Performing Wilcoxon signed-ranks tests across {n_subjects} subjects and {n_nodes} nodes for {track_base}.")

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

    output_dir = os.path.join(root, "group_analysis")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"groupwise_{track_base}_left_vs_right_wilcoxon.csv")
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


import os
import numpy as np
import pandas as pd
from scipy import stats

def groupwise_track_wilcoxon_non_resampled(root, subject_dirs, track_base):
    """
    Performs a Wilcoxon Signed-Ranks Test on the non-resampled (single) mean FA
    values of the left vs. right tract (e.g., AF or CST) across subjects.

    For each subject directory (e.g., /Users/.../Bundles/11),
    it looks for the following files in:
        <subject_dir>/bundle_stats/<track_base>_left_mean_fa.csv
        <subject_dir>/bundle_stats/<track_base>_right_mean_fa.csv

    Each CSV should contain at least one row with a column "Mean_FA".

    The function:
      1. Reads each subject's left/right mean FA,
      2. Performs a paired Wilcoxon test (left vs. right),
      3. Saves a one-row CSV file with the statistic, p-value, number of subjects,
         and group means to:
         <root>/group_analysis/groupwise_<track_base>_left_vs_right_wilcoxon_non_resampled.csv
    """
    left_vals = []
    right_vals = []

    for subj in subject_dirs:
        # Look for the mean FA CSVs in the "bundle_stats" folder
        stats_dir = os.path.join(subj, "bundle_stats")
        left_file = os.path.join(stats_dir, f"{track_base}_left_mean_fa.csv")
        right_file = os.path.join(stats_dir, f"{track_base}_right_mean_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        try:
            df_left = pd.read_csv(left_file)
            df_right = pd.read_csv(right_file)

            # Extract the single mean FA value (assuming column "Mean_FA")
            left_val = df_left["Mean_FA"].iloc[0]
            right_val = df_right["Mean_FA"].iloc[0]

            left_vals.append(left_val)
            right_vals.append(right_val)

        except Exception as e:
            print(f"Error reading files for subject {subj}: {e}. Skipping.")
            continue

    if len(left_vals) == 0:
        print("No valid subjects found for group-wise Wilcoxon test (non-resampled).")
        return

    # Perform the Wilcoxon Signed-Ranks Test
    stat, p_value = stats.wilcoxon(left_vals, right_vals)
    group_mean_left = np.mean(left_vals)
    group_mean_right = np.mean(right_vals)

    # Create the output directory and file path
    output_dir = os.path.join(root, "group_analysis")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(
        output_dir, f"groupwise_{track_base}_left_vs_right_wilcoxon_non_resampled.csv"
    )

    # Prepare a one-row DataFrame with the results
    results_df = pd.DataFrame({
        "Wilcoxon_Statistic": [stat],
        "p_value": [p_value],
        "N_subjects": [len(left_vals)],
        "Group_mean_left": [group_mean_left],
        "Group_mean_right": [group_mean_right]
    })

    # Save results
    results_df.to_csv(output_file, index=False)
    print(f"Group-wise non-resampled Wilcoxon results saved to {output_file}")
    print(f"Wilcoxon statistic: {stat}, p-value: {p_value}")
    print(f"Mean left: {group_mean_left}, Mean right: {group_mean_right}")


def plot_tractometry_results(root, bundles="AF_left AF_right"):
    """
    This functions first creates the necessary TractSeg file with tractometry path, bundles to use
    and a 3D plot option. Also writes the necessary subject info.
    It then creates a png with the txt options.
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


def get_nodewise_fa_range(subject_dirs, track_base="AF"):
    """
    Gathers all node-wise Mean_FA values for the left and right tracts (e.g., AF_left_n100_fa_statistics.csv)
    across all subjects, then computes the global min and max for each side.

    File naming convention per subject:
      <subject_dir>/bundle_stats/<track_base>_left_n100_fa_statistics.csv
      <subject_dir>/bundle_stats/<track_base>_right_n100_fa_statistics.csv

    Each CSV should have columns: ["Node", "Mean_FA", "Std_FA"] (from your run_fa_processing).
    The function:
      - Collects the "Mean_FA" values for each node from all subjects (Left & Right).
      - Returns and prints the min & max for left and right.

    Parameters
    ----------
    subject_dirs : list of str
        Paths to subject directories (e.g., ["/.../Bundles/11", "/.../Bundles/13", ...])
    track_base : str
        The base name of the tract (e.g., "AF" or "CST").

    Returns
    -------
    (left_min, left_max, right_min, right_max) : tuple of floats or None
        Global min and max of Mean_FA across nodes and subjects for Left & Right.
        Returns None if no valid data found.
    """

    all_left_values = []
    all_right_values = []

    for subj in subject_dirs:
        stats_dir = os.path.join(subj, "bundle_stats")
        left_stats_file = os.path.join(stats_dir, f"{track_base}_left_n100_fa_statistics.csv")
        right_stats_file = os.path.join(stats_dir, f"{track_base}_right_n100_fa_statistics.csv")

        # Skip if either file doesn't exist
        if not os.path.exists(left_stats_file) or not os.path.exists(right_stats_file):
            print(f"Missing stats file(s) in {stats_dir}. Skipping subject {subj}.")
            continue

        try:
            df_left = pd.read_csv(left_stats_file)
            df_right = pd.read_csv(right_stats_file)
            # Extend the global lists with all node-wise Mean_FA values
            all_left_values.extend(df_left["Mean_FA"].values)
            all_right_values.extend(df_right["Mean_FA"].values)
        except Exception as e:
            print(f"Error reading stats for subject {subj}: {e}. Skipping.")

    # If no valid data, return None
    if len(all_left_values) == 0 or len(all_right_values) == 0:
        print("No valid node-wise FA data found across subjects.")
        return None

    left_min, left_max = min(all_left_values), max(all_left_values)
    right_min, right_max = min(all_right_values), max(all_right_values)

    # Print results
    print(f"\n--- {track_base} Node-wise FA Range Across All Subjects ---")
    print(f"Left {track_base}:  min = {left_min:.4f}, max = {left_max:.4f}")
    print(f"Right {track_base}: min = {right_min:.4f}, max = {right_max:.4f}\n")

    return left_min, left_max, right_min, right_max


def report_nonresampled_fa_stats(subject_dirs, track_base="AF"):
    """
    Gathers non-resampled mean FA values for left vs. right tracts across subjects,
    computes group mean ± SD, runs a Wilcoxon signed-rank test, and prints a
    summary statement.

    Expects each subject directory to have:
      <subject_dir>/bundle_stats/<track_base>_left_mean_fa.csv
      <subject_dir>/bundle_stats/<track_base>_right_mean_fa.csv

    Each CSV must have a "Mean_FA" column in row 0.

    Parameters
    ----------
    subject_dirs : list of str
        Paths to subject directories (e.g., ['/.../Bundles/11', '/.../Bundles/13', ...])
    track_base : str
        Base name of the tract (e.g., "AF" or "CST").

    Returns
    -------
    None
    """

    left_vals = []
    right_vals = []

    # 1) Collect mean FA from each subject
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

    if len(left_vals) == 0:
        print("No valid subjects found for non-resampled FA. Cannot compute summary.")
        return

    # 2) Compute group mean ± SD
    mean_left = np.mean(left_vals)
    std_left = np.std(left_vals, ddof=1)
    mean_right = np.mean(right_vals)
    std_right = np.std(right_vals, ddof=1)

    # 3) Wilcoxon signed-rank test
    stat, p_value = stats.wilcoxon(left_vals, right_vals)

    # 4) Print the statement
    print(
        f"When averaging along the non-resampled tract, the mean left {track_base} FA was "
        f"{mean_left:.3f} (± {std_left:.3f}) compared with {mean_right:.3f} (± {std_right:.3f}) "
        f"on the right (Wilcoxon signed-rank p={p_value:.3f})."
    )


