import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#from statsmodels.stats.multitest import multipletests
from scipy import stats
from helpers.helpers import run_cmd, create_tractseg_file


def visualize_fa(csv_file):
    """
    Visualizes the FA profile along a tract from a CSV file.
    """

    # Read the FA CSV file into a df
    data = pd.read_csv(csv_file, skiprows=1, header=None, delim_whitespace=True)

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


def process_fa_stats(input_file):
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
    output_dir = os.path.join("root", "group_analysis", "stats")
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


def perform_nodewise_ttest(stats1, stats2, paired=False, fdr_adjust=True):
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

    # Ensure both datasets have the same number of nodes (columns)
    if data1.shape[1] != data2.shape[1]:
        raise ValueError("Both datasets must have the same number of nodes (columns).")

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

    return results


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
    plot_tractometry_results("/media/nas/nikita/test_study2_1sub")

    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/FA/AF_left_n100_fa.csv")
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_fa.csv")
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_right_fa.csv")
    #visualize_peak_length("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_peaks.txt")

    #input_file_left = "/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_left_fa.csv"
    #input_file_right = "/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_right_fa.csv"
    #stats_left = process_fa_stats(input_file_left)
    #stats_right = process_fa_stats(input_file_right)

    #ttest_results = perform_nodewise_ttest(stats_left, stats_right, paired=True, fdr_adjust=True)

    #print(ttest_results.to_string(index=False))
