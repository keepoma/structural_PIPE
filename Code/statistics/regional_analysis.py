import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats
from Code.helpers.helpers import (run_cmd, create_tractseg_file,
                                  get_subject_dirs)


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


def visualize_fa(csv_file):
    """
    Visualizes the FA profile along a tract from a CSV file.
    """

    # Read the FA CSV file into a df
    data = pd.read_csv(csv_file, skiprows=1, header=None, sep='\s+')

    plt.figure(figsize=(10, 5))

    mean_values = data.mean(axis=0)
    std_values = data.std(axis=0)
    x = range(1, len(mean_values) + 1)
    plt.plot(x, mean_values, label="Mean FA", color="blue", linewidth=2)
    plt.fill_between(x, mean_values - std_values, mean_values + std_values,
                     color="blue", alpha=0.3, label="Std. Deviation")

    plt.xlabel("Along-Tract Position (resampled)")
    plt.ylabel("FA Value")
    plt.title(os.path.basename(csv_file))
    plt.legend()
    plt.grid(True)
    plt.show()


def visualize_multi_fa(dirs, csv_files):
    """
    Visualizes multiple FA profiles in one plot
    """
    plt.style.use('fast')
    plt.figure(figsize=(12, 6))

    # Define colors for up to four files
    colors = ["blue", "red", "green", "orange"]

    for idx, csv_file in enumerate(csv_files):
        data = pd.read_csv(csv_file, skiprows=1, header=None, sep='\s+')
        mean_values = data.mean(axis=0)
        std_values = data.std(axis=0)
        x = range(1, len(mean_values) + 1)

        color = colors[idx % len(colors)]
        # Plot mean FA
        plt.plot(x, mean_values, label=f"{os.path.basename(csv_file)} Mean FA", color=color, linewidth=2)
        # Plot standard deviation as a shaded area
        plt.fill_between(x, mean_values - std_values, mean_values + std_values,
                         color=color, alpha=0.3, label=f"{os.path.basename(csv_file)} Std")

    file_names = []
    subject_ids = []
    for file in csv_files:
        file_names.append(os.path.basename(file))
    for dir in dirs:
        subject_ids.append(os.path.basename(dir))

    plt.xlabel("Along-Tract Node Position")
    plt.ylabel("FA Value")
    plt.title(f"FA Profiles for {file_names[0]}, subject {subject_ids[0]} "
              f"and {file_names[1]}, subject {subject_ids[1]}")
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.show()


def visualize_peak_length(peaks_txt):
    """
    Generates mean peak amplitude from 2000 streamlines across 100 points
    """

    # Load data as 1D array
    data_1d = np.loadtxt(peaks_txt)
    n_points = 100
    n_streamlines = 2000

    # Reshape data to (n_streamlines, n_points)
    data_2d = data_1d.reshape(n_streamlines, n_points)

    # Mean and SD plot
    mean_values = data_2d.mean(axis=0)
    std_values = data_2d.std(axis=0)

    plt.figure(figsize=(8, 5))
    plt.plot(mean_values, label="Mean")
    plt.fill_between(
        np.arange(n_points),
        mean_values - std_values,
        mean_values + std_values,
        alpha=0.2, label="±1 SD"
    )

    plt.title("Mean ± SD of Peak Values")
    plt.xlabel("Node index")
    plt.ylabel("Peak amplitude")
    plt.legend()
    plt.show()


def visualize_multi_peak_length(peaks_txt_files):
    """
    Visualizes multiple peak amplitude profiles in one plot.
    """

    plt.style.use('fast')
    plt.figure(figsize=(12, 6))
    colors = ["blue", "red", "green", "orange"]

    n_points = 100
    n_streamlines = 2000

    for idx, txt_file in enumerate(peaks_txt_files):
        data_1d = np.loadtxt(txt_file)
        # Reshape data to (n_streamlines, n_points)
        data_2d = data_1d.reshape(n_streamlines, n_points)

        mean_values = data_2d.mean(axis=0)
        std_values = data_2d.std(axis=0)

        color = colors[idx % len(colors)]
        label_base = os.path.basename(txt_file)

        # Plot the mean
        plt.plot(mean_values, color=color, linewidth=2,
                 label=f"{label_base} Mean")
        # Fill ±1 SD
        x = np.arange(n_points)
        plt.fill_between(x, mean_values - std_values, mean_values + std_values,
                         color=color, alpha=0.2,
                         label=f"{label_base} ±1 SD")

    plt.title("Mean ± SD of Peak Values")
    plt.xlabel("Node Index")
    plt.ylabel("Peak Amplitude")
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.tight_layout()
    plt.show()


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
    visualize_multi_fa(
        subject_dirs,
        AF_n100_files_between_two_subj)

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
