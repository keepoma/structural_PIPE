import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#from statsmodels.stats.multitest import multipletests
from scipy import stats
from contextlib import contextmanager
from helpers import run_cmd, get_subject_paths, get_subject_dirs, get_args

@contextmanager
def change_dir(new_dir):
    """
    Context manager for temporarily changing the working directory
    """

    old_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(old_dir)


def compute_group_response_functions(root, output_dir, nthreads):
    """
    Compute group-average response functions for each tissue type.
    """

    os.makedirs(output_dir, exist_ok=True)
    subject_dirs = get_subject_dirs(root)

    # List of tissue types to process.
    tissue_types = ["wm", "gm", "csf"]
    response_files = {tissue: [] for tissue in tissue_types}

    # Gather response function files for each subject.
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        for tissue in tissue_types:
            tissue_file = os.path.join(paths["five_dwi"], f"{tissue}.txt")
            response_files[tissue].append(tissue_file)

    # Run the responsemean command for each tissue type.
    for tissue in tissue_types:
        group_file = os.path.join(output_dir, f"group_average_response_{tissue}.txt")

        run_cmd([
            "responsemean",
            *response_files[tissue],
            group_file,
            "-nthreads", str(nthreads),
            "-force"
        ])


def population_template_and_register(root, nthreads):
    """
    Build a group-level FOD template from individual subjects’ data,
    then register each subject’s FOD image to that template.
    """

    group_analysis_dir = os.path.join(root, "group_analysis")
    template_dir = os.path.join(group_analysis_dir, "template")
    fod_input_dir = os.path.join(template_dir, "fod_input")
    mask_input_dir = os.path.join(template_dir, "mask_input")
    os.makedirs(fod_input_dir, exist_ok=True)
    os.makedirs(mask_input_dir, exist_ok=True)

    subject_dirs = get_subject_dirs(root)

    # Copy required FOD and mask files to the template input directories.
    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        paths = get_subject_paths(subj_dir)
        files_to_copy = {
            "fod": ("wm_group_average_based_norm.mif", fod_input_dir),
            "mask": ("mask.mif", mask_input_dir)
        }
        for key, (src_filename, dest_dir) in files_to_copy.items():
            src_path = os.path.join(paths["five_dwi"], src_filename)
            dest_path = os.path.join(dest_dir, f"{subject_id}.mif")
            run_cmd(["cp", src_path, dest_path])

    # Build the population FOD template.
    output_template = os.path.join(template_dir, "wmfod_template.mif")
    run_cmd([
        "population_template",
        fod_input_dir,
        "-mask_dir", mask_input_dir,
        output_template,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Register each subject’s FOD image to the template.
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        five_dir = paths["five_dwi"]
        with change_dir(five_dir):
            run_cmd([
                "mrregister",
                "wm_group_average_based_norm.mif",
                "-mask1", "mask.mif",
                os.path.relpath(output_template, five_dir),
                "-nl_warp", "subject2template_warp.mif", "template2subject_warp.mif",
                "-nthreads", str(nthreads),
                "-force"
            ])


def warp_masks_and_create_template_mask(root, nthreads):
    """
    Warp each subject’s mask into template space and compute the intersection
    to form the template mask.
    """

    subject_dirs = get_subject_dirs(root)
    warped_mask_paths = []

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        five_dir = paths["five_dwi"]
        with change_dir(five_dir):
            input_mask = "mask.mif"
            warp_file = "subject2template_warp.mif"
            output_warped_mask = "mask_in_template_space.mif"
            run_cmd([
                "mrtransform", input_mask,
                "-warp", warp_file,
                "-interp", "nearest",
                "-datatype", "bit",
                output_warped_mask,
                "-nthreads", str(nthreads),
                "-force"
            ])
            warped_mask_paths.append(os.path.join(os.getcwd(), output_warped_mask))

    template_mask = os.path.join(root, "group_analysis", "template", "template_mask.mif")

    # Create group-level template mask by taking the min value
    # at each voxel across all subjects masks in template space.
    # It's conservative (a single 0 voxel is enough to exclude)
    run_cmd([
        "mrmath",
        *warped_mask_paths,
        "min",
        template_mask,
        "-datatype", "bit",
        "-nthreads", str(nthreads),
        "-force"
    ])


def create_group_fixel_mask(template_dir, nthreads):
    """
    Create the group fixel mask from the template FOD.
    """

    template_mask = os.path.join(template_dir, "template_mask.mif")
    wm_template = os.path.join(template_dir, "wmfod_template.mif")
    fixel_mask_dir = os.path.join(template_dir, "fixel_mask")
    run_cmd([
        "fod2fixel",
        "-mask", template_mask,
        "-fmls_peak_value", "0.06",
        wm_template,
        fixel_mask_dir, "-force"
    ])


def process_all_subjects_fixels(root):
    """
    Process fixels for each subject:
      - Warp the subject’s FOD into template space
      - Compute fixels from the warped FOD
      - Reorient the fixels
      - Establish fixel correspondence for both FD and FC metrics
      - Compute FC
    """

    subject_dirs = get_subject_dirs(root)
    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        paths = get_subject_paths(subj_dir)
        # Paths in template space.
        template_mask_rel = os.path.join(root, "group_analysis", "template", "template_mask.mif")
        fixel_mask_rel = os.path.join(root, "group_analysis", "template", "fixel_mask")
        fc_dir_rel = os.path.join(root, "group_analysis", "template", "fc")
        # Define a dedicated directory for FD in template space.
        fd_dir = os.path.join(root, "group_analysis", "template", "fd")
        os.makedirs(fd_dir, exist_ok=True)

        with change_dir(paths["five_dwi"]):

            # Warp the FOD into template space.
            run_cmd([
                "mrtransform",
                "wm_group_average_based_norm.mif",
                "-warp", "subject2template_warp.mif",
                "-reorient_fod", "no",
                "fod_in_template_space_not_reoriented.mif",
                "-force"
            ])

            # Compute fixels from the non-reoriented FOD.
            run_cmd([
                "fod2fixel",
                "-mask", template_mask_rel,
                "fod_in_template_space_not_reoriented.mif",
                "fixel_in_template_space_not_reoriented",
                "-afd", "fd.mif",
                "-force"
            ])

            # Reorient the fixels.
            run_cmd([
                "fixelreorient",
                "fixel_in_template_space_not_reoriented",
                "subject2template_warp.mif",
                "fixel_in_template_space",
                "-force"
            ])

            # Establish fixel correspondence for the FD metric.
            run_cmd([
                "fixelcorrespondence",
                os.path.join("fixel_in_template_space", "fd.mif"),
                fixel_mask_rel,
                fd_dir,
                f"{subject_id}.mif",
                "-force"
            ])

            # Compute a fixel metric (FC) using warp2metric.
            run_cmd([
                "warp2metric",
                "subject2template_warp.mif",
                "-fc",
                fixel_mask_rel,
                fc_dir_rel,
                f"{subject_id}.mif",
                "-force"
            ])


def post_process_fixel_metrics(template_dir, subject_dirs):
    """
    Post-process fixel metrics (FC and FDC).
    """

    fc_dir = os.path.join(template_dir, "fc")
    log_fc_dir = os.path.join(template_dir, "log_fc")
    os.makedirs(log_fc_dir, exist_ok=True)

    # Copy required index and directions files.
    run_cmd([
        "cp",
        os.path.join(fc_dir, "index.mif"),
        os.path.join(fc_dir, "directions.mif"),
        log_fc_dir
    ])

    # Compute the logarithm for each subject's fc file.
    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        fc_file = os.path.join(fc_dir, f"{subject_id}.mif")
        log_fc_file = os.path.join(log_fc_dir, f"{subject_id}.mif")
        run_cmd(["mrcalc", fc_file, "-log", log_fc_file, "-force"])

    # Compute the fdc metric.
    fdc_dir = os.path.join(template_dir, "fdc")
    os.makedirs(fdc_dir, exist_ok=True)
    run_cmd(["cp", os.path.join(fc_dir, "index.mif"), fdc_dir])
    run_cmd(["cp", os.path.join(fc_dir, "directions.mif"), fdc_dir])
    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        fd_file = os.path.join(template_dir, "fd", f"{subject_id}.mif")
        fc_file = os.path.join(fc_dir, f"{subject_id}.mif")
        fdc_file = os.path.join(fdc_dir, f"{subject_id}.mif")
        run_cmd(["mrcalc", fd_file, fc_file, "-mult", fdc_file])

def run_group_tractography(template_dir, nthreads):
    """
    Run group-level tractography and subsequent post-processing.
    """

    with change_dir(template_dir):
        # Generate the initial tractography.
        run_cmd([
            "tckgen",
            "-angle", "22.5",
            "-maxlen", "250",
            "-minlen", "10",
            "-power", "1.0",
            "wmfod_template.mif",
            "-seed_image", "template_mask.mif",
            "-mask", "template_mask.mif",
            "-select", "5000000",
            "-cutoff", "0.06",
            "tracks_5_million.tck",
            "-nthreads", str(nthreads),
            "-force"
        ])

        # Apply SIFT to reduce the tractogram.
        run_cmd([
            "tcksift",
            "tracks_5_million.tck",
            "wmfod_template.mif",
            "tracks_2_million_sift.tck",
            "-term_number", "2000000",
            "-nthreads", str(nthreads),
            "-force"
        ])

        # Compute connectivity and perform fixel filtering.
        run_cmd(["fixelconnectivity", "fixel_mask/", "tracks_2_million_sift.tck", "matrix/"])
        run_cmd(["fixelfilter", "fd", "smooth", "fd_smooth", "-matrix", "matrix/"])
        run_cmd(["fixelfilter", "log_fc", "smooth", "log_fc_smooth", "-matrix", "matrix/"])
        run_cmd(["fixelfilter", "fdc", "smooth", "fdc_smooth", "-matrix", "matrix/"])


def create_and_run_glm(template_dir, subject_dirs):
    """
    Create a simple General Lineal Model (GLM) for fixel-based analyses
    and run fixelcfestats for FD, log FC, and FDC.

    This function creates three files:
      - design_matrix.txt: A column of ones (one per subject).
      - contrast_matrix.txt: A single "1" (for a simple effect).
      - files.txt: A list of file names (one per subject), assumed to be named <subject_id>.mif
        in the smoothed metric directories.

    It then runs fixelcfestats for each smoothed metric directory:
      - fd_smooth/
      - log_fc_smooth/
      - fdc_smooth/
    """

    # Define paths for the GLM text files within the template directory.
    design_matrix_path = os.path.join(template_dir, "design_matrix.txt")
    contrast_matrix_path = os.path.join(template_dir, "contrast_matrix.txt")
    files_list_path = os.path.join(template_dir, "files.txt")

    num_subjects = len(subject_dirs)

    # Create the design matrix (each subject gets a "1" on its own line).
    with open(design_matrix_path, "w") as f:
        for _ in range(num_subjects):
            f.write("1\n")

    # Create the contrast matrix (a single "1").
    with open(contrast_matrix_path, "w") as f:
        f.write("1\n")

    # Create the subject files list.
    # Assumes that in each metric's smooth directory, each subject's file is named <subject_id>.mif.
    with open(files_list_path, "w") as f:
        for subj_dir in subject_dirs:
            subject_id = os.path.basename(subj_dir)
            f.write(f"{subject_id}.mif\n")

    # Define the connectivity matrix directory.
    connectivity_matrix = os.path.join(template_dir, "matrix")

    # Run fixelcfestats for each metric.
    # FD
    fd_smooth_dir = os.path.join(template_dir, "fd_smooth")
    stats_fd_dir = os.path.join(template_dir, "stats_fd")
    os.makedirs(stats_fd_dir, exist_ok=True)
    run_cmd([
        "fixelcfestats",
        fd_smooth_dir,
        files_list_path,
        design_matrix_path,
        contrast_matrix_path,
        connectivity_matrix,
        stats_fd_dir,
        "-force"
    ])

    # log FC
    log_fc_smooth_dir = os.path.join(template_dir, "log_fc_smooth")
    stats_log_fc_dir = os.path.join(template_dir, "stats_log_fc")
    os.makedirs(stats_log_fc_dir, exist_ok=True)
    run_cmd([
        "fixelcfestats",
        log_fc_smooth_dir,
        files_list_path,
        design_matrix_path,
        contrast_matrix_path,
        connectivity_matrix,
        stats_log_fc_dir,
        "-force"
    ])

    # FDC
    fdc_smooth_dir = os.path.join(template_dir, "fdc_smooth")
    stats_fdc_dir = os.path.join(template_dir, "stats_fdc")
    os.makedirs(stats_fdc_dir, exist_ok=True)
    run_cmd([
        "fixelcfestats",
        fdc_smooth_dir,
        files_list_path,
        design_matrix_path,
        contrast_matrix_path,
        connectivity_matrix,
        stats_fdc_dir,
        "-force"
    ])


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


def main():
    args = get_args()
    root = os.path.abspath(args.root)
    group_output_directory = os.path.join(root, "group_analysis")
    template_dir = os.path.join(root, "group_analysis", "template")
    subject_dirs = get_subject_dirs(root)
    os.makedirs(group_output_directory, exist_ok=True)

    print(f"\n========= Calculating Group RF =========\n")
    compute_group_response_functions(root, group_output_directory, args.nthreads)

    print(f"\n========= Building FOD template and Registering Subjects =========\n")
    population_template_and_register(root, args.nthreads)

    print(f"\n========= Warping and creating Template Mask =========\n")
    warp_masks_and_create_template_mask(root, args.nthreads)

    print(f"\n========= Creating Group Fixel Mask =========\n")
    create_group_fixel_mask(template_dir, args.nthreads)

    print(f"\n========= Processing fixels for each subject =========\n")
    process_all_subjects_fixels(root)

    print(f"\n========= Post Processing fixel Metrics =========\n")
    post_process_fixel_metrics(template_dir, subject_dirs)

    print(f"\n========= Running group tractography =========\n")
    run_group_tractography(template_dir, args.nthreads)

    print(f"\n========= Create and Run GLM =========\n")
    create_and_run_glm(template_dir, subject_dirs)


if __name__ == "__main__":
    main()
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_fa.csv")
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_right_fa.csv")
    #visualize_peak_length("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_peaks.txt")

    #input_file_left = "/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_left_fa.csv"
    #input_file_right = "/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_right_fa.csv"
    #stats_left = process_fa_stats(input_file_left)
    #stats_right = process_fa_stats(input_file_right)

    #ttest_results = perform_nodewise_ttest(stats_left, stats_right, paired=True, fdr_adjust=True)

    #print(ttest_results.to_string(index=False))
