import os
from contextlib import contextmanager
from first_pipe.helpers import run_cmd, get_subject_paths, get_subject_dirs, get_args

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
    os.makedirs(output_dir, exist_ok=True)
    subject_dirs = get_subject_dirs(root)

    # List of tissue types to process.
    tissue_types = ["wm", "gm", "csf"]
    # Dictionary to hold file paths for each tissue type.
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


def build_and_register_fod_template(root, nthreads, voxel_size="1.75"):
    group_analysis_dir = os.path.join(root, "group_analysis")
    template_dir = os.path.join(group_analysis_dir, "template")
    fod_input_dir = os.path.join(template_dir, "fod_input")
    mask_input_dir = os.path.join(template_dir, "mask_input")

    os.makedirs(fod_input_dir, exist_ok=True)
    os.makedirs(mask_input_dir, exist_ok=True)

    subject_dirs = get_subject_dirs(root)

    # Copy each subject's normalized FOD image and mask into the template directories.
    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        paths = get_subject_paths(subj_dir)

        # Mapping file types to their respective source file names and destination directories.
        files_to_copy = {
            "fod": ("wm_norm.mif", fod_input_dir),
            "mask": ("mask.mif", mask_input_dir)
        }

        for key, (src_filename, dest_dir) in files_to_copy.items():
            src_path = os.path.join(paths["five_dwi"], src_filename)
            dest_path = os.path.join(dest_dir, f"{subject_id}.mif")
            run_cmd(["cp", src_path, dest_path])

    output_template = os.path.join(template_dir, "wmfod_template.mif")
    run_cmd([
        "population_template",
        fod_input_dir,
        "-mask_dir", mask_input_dir,
        output_template,
        "-voxel_size", voxel_size,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Register each subject's FOD image to the template.
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        fod_dir = paths["five_dwi"]
        with change_dir(fod_dir):
            cmd = [
                "mrregister",
                "wm_norm.mif",
                "-mask1", "mask.mif",  # Changed from "dwi_mask_upsampled.mif" to "mask.mif"
                os.path.relpath(output_template, fod_dir),
                "-nl_warp", "subject2template_warp.mif", "template2subject_warp.mif",
                "-nthreads", str(nthreads),
                "-force"
            ]
            run_cmd(cmd)

def warp_masks_and_create_template_mask(root, nthreads):
    subject_dirs = get_subject_dirs(root)
    warped_mask_paths = []

    # Warp each subject's mask into template space.
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        fod_dir = paths["five_dwi"]
        with change_dir(fod_dir):
            input_mask = "mask.mif"  # Changed from "dwi_mask_upsampled.mif" to "mask.mif"
            warp_file = "subject2template_warp.mif"
            output_warped_mask = "mask_in_template_space.mif"
            run_cmd([
                "mrtransform",
                input_mask,
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

if __name__ == "__main__":
    args = get_args()
    root = os.path.abspath(args.root)
    group_output_directory = os.path.join(root, "group_analysis")

    # Compute group-average response functions.
    compute_group_response_functions(root, group_output_directory, args.nthreads)

    # Build the FOD template and register subjects.
    build_and_register_fod_template(root, args.nthreads)

    # Warp masks and compute the intersection template mask.
    warp_masks_and_create_template_mask(root, args.nthreads)
