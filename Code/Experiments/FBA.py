import os
from contextlib import contextmanager
from Code.helpers.helpers import run_cmd, get_subject_paths, get_subject_dirs, get_args
from Code.statistical_analysis import compute_group_response_functions



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