import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from helpers.j_helpers import (run_cmd, create_tractseg_file,
                                  get_subject_dirs)


def process_arcuate_mean_fa(paths, subject_dir):
    """
    Opens the FA CSV files for the left and right arcuate fasciculus,
    computes the mean FA for each tract, and saves the output as two CSV files.

    Expected file paths (within subject_dir):
      - ses_pre/along_tract/FA/AF_left_fa.csv
      - ses_pre/along_tract/FA/AF_right_fa.csv

    Output files (saved in paths["bun_stats_dir"]):
      - AF_left_mean_fa.csv
      - AF_right_mean_fa.csv
    """
    # Define the FA directory for the subject
    fa_dir = os.path.join(subject_dir, "ses_pre", "along_tract", "FA")
    if not os.path.exists(fa_dir):
        print(f"FA directory not found at {fa_dir}.")
        return

    # Ensure the output directory exists
    os.makedirs(paths["bun_stats_dir"], exist_ok=True)

    # List of tracts to process
    tracts = ["AF_left", "AF_right"]

    for tract in tracts:
        # Construct the expected file path for the tract
        file_path = os.path.join(fa_dir, f"{tract}_fa.csv")
        if not os.path.exists(file_path):
            print(f"File {file_path} not found. Skipping {tract}.")
            continue

        # Read FA values from the file
        all_values = []
        with open(file_path, 'r') as f:
            for line in f:
                # Remove null characters and trim whitespace
                line = line.replace('\x00', '').strip()
                if not line or line.startswith('#'):
                    continue
                # Convert values to floats (assuming whitespace-delimited)
                values = [float(x) for x in line.split()]
                all_values.extend(values)

        # Compute the mean FA for this tract
        mean_fa = np.mean(all_values)

        # Create a DataFrame with the result
        df = pd.DataFrame({'Tract': [tract], 'Mean_FA': [mean_fa]})

        # Define output file path
        output_filename = f"{tract}_mean_fa.csv"
        output_file = os.path.join(paths["bun_stats_dir"], output_filename)

        # Save the result as a CSV file
        df.to_csv(output_file, index=False)
        print(f"Mean FA for {tract} saved to {output_file}")


def process_fa_stats(paths, input_file):
    """
    Reads a file containing FA values for streamlines along nodes,
    computes summary statistics, and saves the results into a CSV file under paths["bun_stats_dir"].

    Assumptions:
      - Each row corresponds to a streamline.
      - Each column corresponds to a sampling node.

    The output CSV file will contain one row per node with the following columns:
      - Node (node number, starting at 1)
      - Mean_FA (mean FA value)
      - Std_FA (standard deviation of FA)
    """
    # Read the data from file (ignoring comment lines)
    data = pd.read_csv(input_file, header=None, comment='#', sep=r'\s+')
    data = data.astype(float)

    # Compute node-wise statistics
    mean_per_node = data.mean(axis=0)
    std_per_node = data.std(axis=0)

    # Generate the output filename based on the input file name with .csv extension
    base = os.path.basename(input_file)
    name, _ = os.path.splitext(base)
    output_filename = f"{name}_statistics.csv"

    # Create the output directory
    os.makedirs(paths["bun_stats_dir"], exist_ok=True)
    output_file = os.path.join(paths["bun_stats_dir"], output_filename)

    # Create a DataFrame with the results
    stats_df = pd.DataFrame({
        "Node": mean_per_node.index + 1,
        "Mean_FA": mean_per_node.values,
        "Std_FA": std_per_node.values
    })

    # Save the DataFrame as a CSV file
    stats_df.to_csv(output_file, index=False)
    print(f"Statistics saved to {output_file}")

    # Return a dictionary that includes both the summary stats and the raw data
    return {"mean_per_node": mean_per_node, "std_per_node": std_per_node, "data": data}


def run_fa_processing(paths):
    """
    Processes FA files for specified tracts for the current subject/session.

    Expected FA file names (in the along_tract/FA folder):
      - AF_left_n100_fa.csv
      - AF_right_n100_fa.csv

    This function assumes that 'paths' contains enough information to locate the current session's FA folder.
    Typically, the session directory should be available as a key in 'paths' (e.g., "session_dir" or "atlas_dir").
    """
    # Determine the current session directory.
    # Here, we try to get it from 'session_dir' key; if not available, fall back to 'atlas_dir'.
    session_dir = paths.get("session_dir") or paths.get("atlas_dir")
    if session_dir is None:
        print("Cannot determine the session directory from paths. Skipping FA processing.")
        return

    # Construct the FA directory from the session directory.
    fa_dir = os.path.join(session_dir, "along_tract", "FA")
    if not os.path.exists(fa_dir):
        print(f"FA directory {fa_dir} not found. Skipping processing for current session.")
        return

    tract_names = ["AF_left_n100", "AF_right_n100"]
    for tract in tract_names:
        # Build the expected FA file path for the given tract.
        fa_csv = os.path.join(fa_dir, f"{tract}_fa.csv")
        if os.path.exists(fa_csv):
            print(f"Processing {fa_csv}...")
            # process_fa_stats is assumed to save the stats into the output directory specified in paths.
            process_fa_stats(paths, fa_csv)
        else:
            print(f"File {fa_csv} not found in {fa_dir}. Skipping.")


def groupwise_af_ttest(root, subject_dirs):
    """
    Performs a node-wise paired t-test comparing FA values for the left and right
    arcuate fasciculus across subjects (using only the ses_pre session).

    For each subject, it expects two whitespace-delimited files:
      - <subject_dir>/ses_pre/along_tract/FA/AF_left_n100_fa.csv
      - <subject_dir>/ses_pre/along_tract/FA/AF_right_n100_fa.csv

    Each file should have 2000 rows (streamlines) and 100 columns (nodes).

    The function computes the mean FA per node (averaging over streamlines) for each file,
    performs a node-wise paired t-test (comparing left vs. right across subjects),
    applies FDR correction, and saves the results as:
        <root>/group_analysis/groupwise_AF_left_vs_right_ttest.csv
    """
    left_list = []
    right_list = []

    for subj in subject_dirs:
        # Only use the ses_pre session
        fa_dir = os.path.join(subj, "ses_pre", "along_tract", "FA")
        left_file = os.path.join(fa_dir, "AF_left_n100_fa.csv")
        right_file = os.path.join(fa_dir, "AF_right_n100_fa.csv")

        if not os.path.exists(left_file) or not os.path.exists(right_file):
            print(f"Missing file(s) for subject {subj}. Skipping.")
            continue

        try:
            # Use delim_whitespace=True so that the 2000×100 structure is properly parsed.
            left_df = pd.read_csv(left_file, comment='#', header=None, sep='\s+')
            right_df = pd.read_csv(right_file, comment='#', header=None, sep='\s+')
        except Exception as e:
            print(f"Error reading files for subject {subj}: {e}. Skipping.")
            continue

        # Check that each file has 100 columns (nodes)
        if left_df.shape[1] != 100 or right_df.shape[1] != 100:
            print(f"Subject {subj} does not have 100 nodes. Skipping. "
                  f"Left shape: {left_df.shape}, Right shape: {right_df.shape}")
            continue

        # Compute the mean FA per node (average over the 2000 streamlines = rows)
        left_mean = left_df.mean(axis=0).values   # Shape: (100,)
        right_mean = right_df.mean(axis=0).values

        left_list.append(left_mean)
        right_list.append(right_mean)

    if len(left_list) == 0:
        print("No valid subjects found for group-wise t-test.")
        return

    # Convert the lists to arrays of shape (n_subjects, 100)
    left_array = np.array(left_list)
    right_array = np.array(right_list)
    n_subjects, n_nodes = left_array.shape
    print(f"Performing paired t-tests across {n_subjects} subjects and {n_nodes} nodes.")

    t_stats = []
    p_values = []
    for node in range(n_nodes):
        t, p = stats.ttest_rel(left_array[:, node], right_array[:, node])
        t_stats.append(t)
        p_values.append(p)

    # FDR correction using Benjamini-Hochberg
    _, p_values_fdr, _, _ = multipletests(p_values, method='fdr_bh')

    results_df = pd.DataFrame({
        "Node": np.arange(1, n_nodes + 1),
        "t_stat": t_stats,
        "p_value": p_values,
        "p_value_fdr": p_values_fdr
    })

    output_dir = os.path.join(root, "group_analysis")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "groupwise_AF_left_vs_right_ttest.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Group-wise paired t-test results saved to {output_file}")


def plot_tractometry_results(root, bundles="AF_left AF_right CC_5 CC_6 SCP_left"):
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


if __name__ == "__main__":
    root = "/media/nas/nikita/test_study2_1sub"
    subject_dirs = get_subject_dirs(root)
    #plot_tractometry_results("/media/nas/nikita/test_study2_1sub")

    AF_n100_peaks = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/peaks/AF_left_n100_peaks.txt",
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/peaks/AF_right_n100_peaks.txt"
    )
    #visualize_multi_peak_length(AF_peaks)

    AF_n100_fa_values = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_n100_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_right_n100_fa.csv"
    )
    AF_n100_files_between_two_subj = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_n100_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_588/along_tract/FA/AF_left_n100_fa.csv"
    )
    #visualize_multi_fa(AF_fa_values)

    """
    visualize_multi_fa(
        subject_dirs,
        AF_n100_files_between_two_subj)
    
    input_file_left = AF_n100_fa_values[0]
    input_file_right = AF_n100_fa_values[1]
    stats_n100_AF_left = process_fa_stats(root, input_file_left)
    stats_n100_AF_right = process_fa_stats(root, input_file_right)
    ttest_results = perform_nodewise_ttest(
        stats_n100_AF_left,
        stats_n100_AF_right,
        file1=input_file_left,
        file2=input_file_right,
        paired=True,
        fdr_adjust=True)
    """


    AF_fa_values = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_right_fa.csv"
    )
