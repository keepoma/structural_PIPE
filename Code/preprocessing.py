#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse


def run_cmd(cmd, verbose=True):
    """
    Helper function to run a system command via subprocess.
    May print the command and is able to raise errors.
    """
    if verbose:
        print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def process_subject(subject_dir, nthreads):
    """
    Runs the entire preproc pipeline for one subject.

    Parameters:
    -----------
    subject_dir: str
        Path to the subject directory, e.g. /path/to/study/subject_001
    nthreads: int
        Number of threads to pass to MRtrix commands (via -nthreads). Will try to use min 4
    """
    # Expect the raw data in /path/to/study/subject_001/raw/1_raw/
    one_raw = os.path.join(subject_dir, "raw", "1_raw")
    if not os.path.isdir(one_raw):
        print(f"Warning: {one_raw} does not exist! Skipping.")
        return

    # Create 2_nifti and 5_dwi directories inside raw
    two_nifti = os.path.join(subject_dir, "raw", "2_nifti")
    five_dwi = os.path.join(subject_dir, "raw", "5_dwi")

    os.makedirs(two_nifti, exist_ok=True)
    os.makedirs(five_dwi, exist_ok=True)

    # Helper function for building paths
    def one_path(subpath):
        return os.path.join(one_raw, subpath)

    # 1. Convert scans to NIFTI or MIF
    # Adjust if your input file names differ
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path("006_T1w_MPR"), os.path.join(two_nifti, "t1.nii.gz"),
        "-force"
    ])
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path("008_T2w_SPC"), os.path.join(two_nifti, "t2.nii.gz"),
        "-force"
    ])
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path("009_t2_space_dark-fluid_sag_p2_iso_0_8"),
        os.path.join(two_nifti, "t2_df.nii.gz"),
        "-force"
    ])
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        one_path("016_dMRI_dir98_AP"), os.path.join(five_dwi, "dwi_ap.mif"),
        "-force"
    ])
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        one_path("019_dMRI_dir98_PA"), os.path.join(five_dwi, "dwi_pa.mif"),
        "-force"
    ])

    # 2. Preprocessing steps
    def five_path(subpath):
        return os.path.join(five_dwi, subpath)

    # Combine AP/PA
    run_cmd([
        "mrcat", "-nthreads", str(nthreads),
        five_path("dwi_ap.mif"), five_path("dwi_pa.mif"), five_path("dwi_all.mif"),
        "-axis", "3"
    ])

    # Denoise
    run_cmd([
        "dwidenoise", "-nthreads", str(nthreads),
        five_path("dwi_all.mif"), five_path("dwi_den.mif"),
        "-noise", five_path("noise.mif")
    ])

    # Residual
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        five_path("dwi_all.mif"), five_path("dwi_den.mif"),
        "-subtract", five_path("residual.mif")
    ])

    # Gibbs unringing
    run_cmd([
        "mrdegibbs", "-nthreads", str(nthreads),
        five_path("dwi_den.mif"), five_path("dwi_den_unr.mif"),
        "-axes", "0,1"
    ])

    # Motion/distortion correction with FSL
    run_cmd([
        "dwifslpreproc", "-nthreads", str(nthreads),
        "-rpe_all", "-pe_dir", "AP",
        five_path("dwi_den_unr.mif"), five_path("dwi_den_unr_pre.mif")
    ])

    # Bias field correction
    run_cmd([
        "dwibiascorrect", "-nthreads", str(nthreads),
        "ants",
        five_path("dwi_den_unr_pre.mif"), five_path("dwi_den_unr_pre_unbia.mif"),
        "-bias", five_path("bias.mif")
    ])

    # Create mask
    run_cmd([
        "dwi2mask", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif")
    ])

    # Skull stripping
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif"),
        "-mult", five_path("dwi_den_unr_pre_unbia_skull.mif")
    ])

    # Response estimation
    run_cmd([
        "dwi2response", "-nthreads", str(nthreads),
        "dhollander",
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("gm.txt"), five_path("csf.txt")
    ])

    # FOD estimation
    run_cmd([
        "dwi2fod", "-nthreads", str(nthreads),
        "msmt_csd",
        "-mask", five_path("mask.mif"),
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("wm.mif"),
        five_path("gm.txt"), five_path("gm.mif"),
        five_path("csf.txt"), five_path("csf.mif")
    ])

    # Intensity normalization
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        five_path("wm.mif"), five_path("wm_norm.mif"),
        five_path("gm.mif"), five_path("gm_norm.mif"),
        five_path("csf.mif"), five_path("csf_norm.mif"),
        "-mask", five_path("mask.mif")
    ])


def main():
    parser = argparse.ArgumentParser(
        description="Run MRtrix pipeline for all subjects in a root directory."
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Path to the root folder containing subject_xxx subdirectories."
    )
    parser.add_argument(
        "--nthreads",
        type=int,
        default=max(4, os.cpu_count() - 10),
        help="Number of threads to pass to MRtrix commands."
    )
    args = parser.parse_args()

    # Collect all subject directories under root
    # e.g., root/subject_001, root/subject_002, etc.
    root = os.path.abspath(args.root)
    subject_dirs = sorted([
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    ])

    # Run the pipeline for each subject
    for subj_dir in subject_dirs:
        print(f"\n=== Processing subject: {os.path.basename(subj_dir)} ===")
        process_subject(subj_dir, args.nthreads)


if __name__ == "__main__":
    main()
