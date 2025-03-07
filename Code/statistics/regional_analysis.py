import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats
from helpers.helpers import (run_cmd, create_tractseg_file,
                                  get_subject_dirs)
from regional_visualization import visualize_multi_peak_length_side_by_side


def compute_mean_fa(csv_file):
    """
    Reads a FA file.
    Coded for non-resampled streamlines, so unequal length is expected (non-matchign columns)
    """

    all_values = []
    with open(csv_file, 'r') as f:
        for line in f:

            # I ran into some \x00 characters within some files, so im removing them
            line = line.replace('\x00', '')

            # Strip whitespace and skip comment
            line = line.strip()
            if line.startswith('#'):
                continue
            # Split on whitespace, convert to floats
            values = [float(x) for x in line.split()]
            # Extend the main list
            all_values.extend(values)

    # Compute mean over the entire set of values
    mean_fa = np.mean(all_values)
    return mean_fa


def process_fa_stats(root, input_file):
    """
    Reads a file containing FA values for streamlines along nodes,
    computes summary statistics, and saves the results into a text file under
    root/group_analysis/stats.

    Assumptions:
      - Each row corresponds to a streamline.
      - Each column corresponds to a sampling node.

    The output file will contain one row per node with the following columns:
      - Node (node number, starting at 1)
      - Mean_FA (mean FA value)
      - Std_FA (standard deviation of FA)
    """

    data = pd.read_csv(input_file, header=None, comment='#', sep=r'\s+')
    data = data.astype(float)

    # Compute node-wise mean and standard deviation
    mean_per_node = data.mean(axis=0)
    std_per_node = data.std(axis=0)

    # Generate the output filename based on the input file name
    base = os.path.basename(input_file)
    name, _ = os.path.splitext(base)
    output_filename = f"{name}_statistics.txt"

    # Define the target directory
    tractseg_tractometry_dir = os.path.join(root, "group_analysis", "tractseg_tractometry")
    os.makedirs(tractseg_tractometry_dir, exist_ok=True)
    output_dir = os.path.join(tractseg_tractometry_dir, "FA_stats")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, output_filename)

    # Create a DF with the results
    stats_df = pd.DataFrame({
        "Node": mean_per_node.index + 1,
        "Mean_FA": mean_per_node.values,
        "Std_FA": std_per_node.values
    })

    # Save the DataFrame to a text file with tab-separated columns.
    stats_df.to_csv(output_file, sep='\t', index=False)
    print(f"Statistics saved to {output_file}")

    # Return a dictionary that includes both the summary stats and the raw data
    return {"mean_per_node": mean_per_node, "std_per_node": std_per_node, "data": data}


def perform_nodewise_ttest(stats1, stats2, file1, file2, paired=False, fdr_adjust=True):
    """
    Performs a t-test at each node between two FA datasets
    and applies FDR correction to the p-values.

    Parameters:
        stats1, stats2 (dict): Outputs from process_fa_stats containing the "data" field.
        paired (bool): If True, a paired t-test is conducted (ttest_rel); otherwise, an independent t-test (ttest_ind) is used.
        fdr_adjust (bool): If True, FDR correction is applied to the p-values.

    Returns:
        results (DataFrame): Contains for each node:
            - Node: Node index
            - t_stat: The t statistic
            - p_value: The original p-value from the t-test
            - p_value_fdr: The FDR-adjusted p-value (if fdr_adjust is True; otherwise same as p_value)
    """

    # Extract the data matrices from the stats dictionaries.
    data1 = stats1["data"]
    data2 = stats2["data"]
    t_stats = []
    p_values = []
    nodes = range(1, data1.shape[1] + 1)

    # Conduct the t-test at each node (adjust index by subtracting 1 for actual column position)
    for i, node in enumerate(nodes):
        values1 = data1.iloc[:, i]
        values2 = data2.iloc[:, i]
        if paired:
            t_stat, p_val = stats.ttest_rel(values1, values2)
        else:
            t_stat, p_val = stats.ttest_ind(values1, values2)
        t_stats.append(t_stat)
        p_values.append(p_val)

    # Apply FDR correction using the Benjamini-Hochberg method if requested
    if fdr_adjust:
        # multipletests returns a tuple where the adjusted p-values are at index 1
        _, p_values_fdr, _, _ = multipletests(p_values, method='fdr_bh')
    else:
        p_values_fdr = p_values

    results = pd.DataFrame({
        "Node": list(nodes),
        "t_stat": t_stats,
        "p_value": p_values,
        "p_value_fdr": p_values_fdr
    })

    tractseg_tractometry_dir = os.path.join(root, "group_analysis", "tractseg_tractometry")
    output_dir = os.path.join(tractseg_tractometry_dir, "FA_stats")

    # Create a filename
    basename1 = os.path.splitext(os.path.basename(file1))[0]
    basename2 = os.path.splitext(os.path.basename(file2))[0]
    paired_str = "paired" if paired else "unpaired"
    fdr_str = "with_FDR" if fdr_adjust else "no_FDR"

    filename = f"t_test_{basename1}_{basename2}_{paired_str}_{fdr_str}.txt"
    filepath = os.path.join(output_dir, filename)

    # Save to a tab-delimited text file
    results.to_csv(filepath, sep=",", index=False)


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
    #root = "/media/nas/nikita/test_study2_1sub"
    #subject_dirs = get_subject_dirs(root)
    #plot_tractometry_results("/media/nas/nikita/test_study2_1sub")

    AF_n100_peaks = (
        "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/along_tract/AF_left_peaks.txt",
        "/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me/along_tract/AF_right_peaks.txt"
    )
    visualize_multi_peak_length_side_by_side(AF_n100_peaks)

    """
        AF_n100_fa_values = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_n100_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_right_n100_fa.csv"
    )
    AF_n100_files_between_two_subj = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_n100_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_588/along_tract/FA/AF_left_n100_fa.csv"
    )
    visualize_multi_fa(AF_fa_values)
    visualize_multi_fa(subject_dirs, AF_n100_files_between_two_subj)
    """


    """
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

    """
        AF_fa_values = (
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_fa.csv",
        "/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_right_fa.csv"
    )
    mean_r_left_af = compute_mean_fa(AF_n100_fa_values[0])
    mean_r_right_af = compute_mean_fa(AF_n100_fa_values[1])
    mean_left_af = compute_mean_fa(AF_fa_values[0])
    mean_right_af = compute_mean_fa(AF_fa_values[1])
    print("mean_r_left_af", mean_r_left_af)
    print("mean_r_right_af", mean_r_right_af)
    print("mean_left_af", mean_left_af)
    print("mean_right_af", mean_right_af)
    """

