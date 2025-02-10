import os
from first_pipe.helpers import run_cmd, get_subject_paths, get_args

"""
Code intended for statistical analysis
1. Conduct fixel-based analysis
"""

# Estimate mean response function across whole group
def compute_group_response_functions(root, output_dir, nthreads):

    os.makedirs(output_dir, exist_ok=True)

    # Get a sorted list of subject directories under the given root and
    # Exclude any folder that should not be considered a subject (e.g., "group_analysis")
    subject_dirs = sorted([
        os.path.join(root, d)
        for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)) and d != "group_analysis"
    ])

    # Lists to collect individual response function file paths.
    wm_files = []
    gm_files = []
    csf_files = []

    # Iterate over each subject and retrieve response function files.
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)

        wm_file = os.path.join(paths["five_dwi"], "wm.txt")
        gm_file = os.path.join(paths["five_dwi"], "gm.txt")
        csf_file = os.path.join(paths["five_dwi"], "csf.txt")

        if os.path.isfile(wm_file):
            wm_files.append(wm_file)
        else:
            print(f"Warning: {wm_file} not found in subject {subj_dir}")

        if os.path.isfile(gm_file):
            gm_files.append(gm_file)
        else:
            print(f"Warning: {gm_file} not found in subject {subj_dir}")

        if os.path.isfile(csf_file):
            csf_files.append(csf_file)
        else:
            print(f"Warning: {csf_file} not found in subject {subj_dir}")

    # Define output file paths for the group-average response functions.
    group_wm = os.path.join(output_dir, "group_average_response_wm.txt")
    group_gm = os.path.join(output_dir, "group_average_response_gm.txt")
    group_csf = os.path.join(output_dir, "group_average_response_csf.txt")

    # Build and run the MRtrix3 'responsemean' command for WM.
    if wm_files:
        run_cmd([
            "responsemean",
            *wm_files,
            group_wm,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print("No WM response files found.")

    # Build and run the MRtrix3 'responsemean' command for GM.
    if gm_files:
        run_cmd([
            "responsemean",
            *gm_files,
            group_gm,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print("No GM response files found.")

    # Build and run the MRtrix3 'responsemean' command for CSF.
    if csf_files:
        run_cmd([
            "responsemean",
            *csf_files,
            group_csf,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print("No CSF response files found.")

    print("Group-average response functions computed and saved.")


if __name__ == "__main__":
    args = get_args()
    root = os.path.abspath(args.root)
    group_output_directory = os.path.join(root, "group_analysis")
    compute_group_response_functions(root, group_output_directory, args.nthreads)
