import glob
import os
import argparse
import pandas as pd
from helpers import run_cmd, get_args

"""
This code is also able to run as a standalone if the preprocessing part is
to be skipped
"""

def process_subject(subject_path, tract_names, nthreads=max(4, os.cpu_count() - 10)):
    """
    Process a single subject by looping through each tract.
    This is the only function imported from this module into main_first_pipe
    """

    for tract_name in tract_names:

        # Build file paths for segmentation files and outputs
        bundle_path = os.path.join(subject_path, "tractseg_output", "bundle_segmentations",
                                   f"{tract_name}.nii.gz")
        bundle_b_path = os.path.join(subject_path, "tractseg_output", "endings_segmentations",
                                     f"{tract_name}_b.nii.gz")
        bundle_e_path = os.path.join(subject_path, "tractseg_output", "endings_segmentations",
                                     f"{tract_name}_e.nii.gz")

        # Output tracking paths
        tracking_dir = os.path.join(subject_path, "tractseg_output", "FOD_iFOD2_trackings")
        os.makedirs(tracking_dir, exist_ok=True)
        tck_path = os.path.join(tracking_dir, f"{tract_name}.tck")
        tck_N100 = os.path.join(tracking_dir, f"{tract_name}_N100.tck")

        # Along-tract output directory
        along_dir = os.path.join(subject_path, "along_tract")
        os.makedirs(along_dir, exist_ok=True)
        adc_csv = os.path.join(along_dir, f"{tract_name}_adc.csv")
        fa_csv = os.path.join(along_dir, f"{tract_name}_fa.csv")
        #c_csv = os.path.join(along_dir, f"{tract_name}_c.csv")

        # 1. Track generation using tckgen
        run_cmd([
            "tckgen",
            "-algorithm", "iFOD2",
            os.path.join(subject_path, "raw", "5_dwi", "wm_norm.mif"),
            tck_path,
            "-seed_image", bundle_path,
            "-mask", bundle_path,
            "-include", bundle_b_path,
            "-include", bundle_e_path,
            "-minlength", "40",
            "-maxlength", "250",
            "-seeds", "1000000",
            "-select", "2000",
            "-cutoff", "0.05",
            "-nthreads", str(nthreads),
            "-force"
        ])

        # 2. Resampling using tckresample
        run_cmd([
            "tckresample",
            tck_path,
            "-num_points", "100",
            tck_N100,
            "-nthreads", str(nthreads),
            "-force"
        ])

        # 3. Sample ADC values along the tract
        run_cmd([
            "tcksample",
            tck_N100,
            os.path.join(subject_path, "raw", "5_dwi", "adc.mif"),
            adc_csv,
            "-fo"
        ])

        # 4. Sample FA values along the tract
        run_cmd([
            "tcksample",
            tck_N100,
            os.path.join(subject_path, "raw", "5_dwi", "fa.mif"),
            fa_csv,
            "-fo"
        ])


def process_all_subjects(root, tract_names_file, nthreads=max(4, os.cpu_count() - 10)):
    """
    Process all subject directories under the root directory.
    """

    # Load tract names from the CSV file; expecting one tract name per row (without header)
    tract_names = pd.read_csv(tract_names_file, header=None)[0].tolist()

    # Loop over each subject (only directories) in the given root
    subject_paths = glob.glob(os.path.join(root, '*'))
    for subject_path in subject_paths:
        if os.path.isdir(subject_path):
            print(f"\n===== Applying segmentation: {subject_path} =====")
            process_subject(subject_path, tract_names, nthreads)


def main():
    args = get_args()
    process_all_subjects(args.root, "tract_name.txt", args.nthreads)


if __name__ == "__main__":
    main()
