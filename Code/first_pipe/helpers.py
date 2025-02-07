import os
import subprocess


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
        "mat_dir": os.path.join(subject_dir, "mat")
    }
    return paths
