import os
import subprocess
import argparse


def run_cmd(cmd):
    """
    Runs a system command via subprocess.
    Prints the command and raises an error if the command fails.
    """
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def get_subject_paths(subject_dir):
    """
    Returns a dictionary of standardized paths for the given subject directory.
    """

    paths = {
        "raw_dir": os.path.join(subject_dir, "raw"),
        "one_raw": os.path.join(subject_dir, "raw", "1_raw"),
        "two_nifti": os.path.join(subject_dir, "raw", "2_nifti"),
        "five_dwi": os.path.join(subject_dir, "raw", "5_dwi"),
        "mat_dir": os.path.join(subject_dir, "mat"),
        "tck_dir": os.path.join(subject_dir, "tck"),
        "atlas_dir": os.path.join(subject_dir, "atlas")
    }
    return paths


def get_subject_dirs(root, exclude="group_analysis"):
    """
    Return a sorted list of subject directories excluding a specific folder
    """

    subject_dirs = sorted([
        os.path.join(root, d)
        for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)) and d != exclude
    ])
    return subject_dirs


def get_args():
    """
    Sets up and returns the parsed command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Run pipeline for all subjects in a root directory. WILL OVERWRITE FILES"
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Path to the root folder containing subject subdirectories."
    )
    parser.add_argument(
        "--nthreads",
        type=int,
        default=max(4, os.cpu_count() - 10),
        help=("Number of threads to pass to MRtrix commands. Will attempt to use "
              "max available threads - 10, if not possible attempts 4.")
    )
    return parser.parse_args()


def prompt_for_folder(default, description):
    """
    Prompts the user for a folder name with a given description.
    Returns the user's input or the default if nothing is entered.
    """

    user_input = input(f"Enter folder name for {description} [default: {default}] (press Enter for default): ").strip()
    return user_input if user_input else default