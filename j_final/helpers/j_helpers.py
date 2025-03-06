import os
import subprocess
import argparse
import logging
from datetime import datetime


def run_cmd(cmd):
    """
    Runs a system command via subprocess.
    Prints the command and raises an error if the command fails.
    """
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def get_subject_paths(subject_dir, session):
    """
    Returns a dictionary of standardized paths for the given subject directory.
    """

    session_path = os.path.join(subject_dir, session)
    paths = {
        "anat_dir": os.path.join(session_path, "anat"),
        "dwi_dir": os.path.join(session_path, "dwi"),
        "fmap_dir": os.path.join(session_path, "fmap"),
        "func_dir": os.path.join(session_path, "func"),
        "mat_dir": os.path.join(subject_dir, "mat"),
        "tck_dir": os.path.join(subject_dir, "tck"),
        "atlas_dir": os.path.join(subject_dir, "atlas"),
        "connectome_dir": os.path.join(subject_dir, "connectome"),
        "tractseg_dir": os.path.join(subject_dir, "tractseg_output")
    }
    return paths


def get_subject_dirs(root):
    """
    Return a sorted list of subject directories excluding a specific folder
    """

    subject_dirs = sorted([
        os.path.join(root, d)
        for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)) and d not in ("group_analysis", "logs")
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
        default=max(4, os.cpu_count() - 50),
        help=("Number of threads to pass to MRtrix commands. Will attempt to use "
              "max available threads - 50, if not possible attempts 4.")
    )
    return parser.parse_args()


def prompt_for_folder(default, description):
    """
    Prompts the user for a folder name with a given description.
    Returns the user's input or the default if nothing is entered.
    """

    user_input = input(f"Enter folder name for {description} [default: {default}] (press Enter for default): ").strip()
    return user_input if user_input else default


def ask_yes_no(question):
    """
    Prompt the user with a yes/no question until
    they enter a valid response. Returns True or False
    """

    while True:
        answer = input(question + " [y/n]: ").strip().lower()
        if answer in ("y", "yes"):
            return True
        elif answer in ("n", "no"):
            return False
        else:
            print("Invalid input. Please type 'y' or 'n'.")

def fancy_print(action, subj_dir):
    """
    Prints a standardized header
    """

    subject_name = os.path.basename(subj_dir)
    print(f"\n========= {action} for Subject: {subject_name} =========\n")


def logs(root):
    # Create the logs folder under args.root
    log_dir = os.path.join(root, "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Create a timestamped log filename
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = os.path.join(log_dir, f"log_{timestamp}.log")

    # Configure logging to write to the log file
    logging.basicConfig(
        filename=log_filename,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s'
    )


def create_tractseg_file(root, tractometry_path, bundles, plot3d, output_file, group=0, age=0.0, sex=0):
    """
    Create a subjects file for plot_tractometry_results.

    The second column has to be 'group' (for a group comparison; containing only 0 or 1) or
    'target' (for a correlation analysis; containing the value you want to calculate the correlation for).
    """

    # Header lines for the subjects file.
    header_lines = [
        f"# tractometry_path={tractometry_path}",
        f"# bundles={bundles}",
        f"# plot_3D={plot3d}",
        "",
        "subject_id group Age Sex"
    ]

    # Retrieve subject IDs from directories under subjects_root.
    subject_dirs = get_subject_dirs(root)

    # Create a line for each subject with default values (group=0, Age=0.0, Sex=0).
    subject_lines = [f"{os.path.basename(subject)} {group} {age} {sex}" for subject in subject_dirs]

    # Write the header and subject lines to the output file.
    with open(output_file, "w") as f:
        f.write("\n".join(header_lines) + "\n")
        f.write("\n".join(subject_lines) + "\n")

    print(f"Subjects file created: {output_file}")
    return output_file


if __name__ == '__main__':
    create_tractseg_file(
        root="/media/nas/nikita/test_study2_1sub",
        tractometry_path="/media/nas/nikita/test_study2_1sub/group_analysis/tractseg_tractometry/SUBJECT_ID/Tractometry.csv",
        bundles="AF_left AF_right CC_5 CC_6 SCP_left",
        plot3d="/media/nas/nikita/test_study2_1sub/test_302/tractseg_output",
        output_file="/media/nas/nikita/test_study2_1sub/group_analysis/tractseg_tractometry/subjects.txt"
    )


