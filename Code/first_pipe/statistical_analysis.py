import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_rel, ttest_ind
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
        fod_dir = paths["five_dwi"]
        with change_dir(fod_dir):
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
      - Compute a fixel metric (FC) using warp2metric.
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

            # Establish fixel correspondence for the FC metric.
            run_cmd([
                "fixelcorrespondence",
                os.path.join("fixel_in_template_space", "fd.mif"),
                fixel_mask_rel,
                fc_dir_rel,
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
                "IN.mif",
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
            "tracks_20_million.tck",
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


def calculate_fa_stats(csv_file):
    """
    Reads a CSV file containing FA values along a tract and computes summary statistics.

    Assumes that each row corresponds to one subject and each column corresponds to an
    along-tract position.

    Returns:
        stats_dict (dict): Dictionary containing:
            - "mean_per_position": Series with the mean FA for each along-tract position.
            - "std_per_position": Series with the standard deviation of FA per position.
            - "overall_mean": Overall mean FA across all positions.
            - "overall_std": Overall standard deviation of FA values.
            - "data": The original DataFrame.
    """
    # Read the CSV file into a DataFrame
    data = pd.read_csv(csv_file, skiprows=1, header=None, delim_whitespace=True)
    data = data.astype(float)

    # Calculate mean and standard deviation for each along-tract position (i.e., each column)
    mean_per_position = data.mean(axis=0)
    std_per_position = data.std(axis=0)

    # Calculate overall mean and standard deviation across all FA values in the DataFrame
    overall_mean = data.values.mean()
    overall_std = data.values.std()

    stats_dict = {
        "mean_per_position": mean_per_position,
        "std_per_position": std_per_position,
        "overall_mean": overall_mean,
        "overall_std": overall_std,
        "data": data
    }

    return stats_dict


def perform_ttest(stats1, stats2, paired=True):
    """
    Performs a t-test comparing the per-subject mean FA values from two FA statistics dictionaries.
    """

    # Compute per-subject mean FA values from the original data
    subject_means1 = stats1["data"].mean(axis=1)
    subject_means2 = stats2["data"].mean(axis=1)

    # For a paired t-test, ensure both arrays have the same length.
    if paired and len(subject_means1) != len(subject_means2):
        raise ValueError("For a paired t-test, both datasets must have the same number of subjects.")

    if paired:
        t_stat, p_value = ttest_rel(subject_means1, subject_means2)
    else:
        t_stat, p_value = ttest_ind(subject_means1, subject_means2)

    return {"t_statistic": t_stat, "p_value": p_value}


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
    #compute_group_response_functions(root, group_output_directory, args.nthreads)

    print(f"\n========= Building FOD template and Registering Subjects =========\n")
    #population_template_and_register(root, args.nthreads)

    print(f"\n========= Warping and creating Template Mask =========\n")
    # warp_masks_and_create_template_mask(root, args.nthreads)

    print(f"\n========= Creating Group Fixel Mask =========\n")
    #create_group_fixel_mask(template_dir, args.nthreads)

    print(f"\n========= Processing fixels for each subject =========\n")
    #process_all_subjects_fixels(root)

    print(f"\n========= Post Processing fixel Metrics =========\n")
    #post_process_fixel_metrics(template_dir, subject_dirs)

    print(f"\n========= Running group tractography =========\n")
    run_group_tractography(template_dir, args.nthreads)


if __name__ == "__main__":
    main()
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_fa.csv")
    #visualize_fa("/media/nas/nikita/test_study2_1sub/test_302/along_tract/AF_right_fa.csv")
    #visualize_peak_length("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_peaks.txt")
    #print(calculate_fa_stats("/media/nas/nikita/test_study2_1sub/test_302/along_tract/CST_left_peaks.txt"))
