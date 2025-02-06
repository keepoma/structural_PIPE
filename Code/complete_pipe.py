import os
import subprocess
import argparse
from bundle_pipe import process_subject
import pandas as pd

"""
General overview and reasoning behind this code:
1) Switched from os.system() to subprocess for better output/error handling
2) Implemented argparse for command-line arguments to make code more versatile
3) Defined functions for code readability and reusability
4) Import process_subject function from bundle_pipe.py. Keeps the code clean and
facilitates troubleshooting
"""

def run_cmd(cmd):
    """
    Function to run a system command via subprocess.
    Prints the command and is able to raise errors.
    """
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def preprocess_subject(subject_dir, nthreads):
    """
    Runs the entire preproc pipeline for one subject.
    subject_dir: str, expects path to the subject directory and
    raw data in /path/to/study/subject_001/raw/1_raw/
    nthreads: int, number of threads to pass to MRtrix commands
    """
    # Simplify path to 1_raw and check that it exists, stops execution if path not present
    one_raw = os.path.join(subject_dir, "raw", "1_raw")
    if not os.path.isdir(one_raw):
        print(f"Warning: {one_raw} does not exist! Exiting.")
        raise SystemExit(1)

    # Create 2_nifti and 5_dwi directories inside raw
    two_nifti = os.path.join(subject_dir, "raw", "2_nifti")
    five_dwi = os.path.join(subject_dir, "raw", "5_dwi")
    os.makedirs(two_nifti, exist_ok=True)
    os.makedirs(five_dwi, exist_ok=True)

    # FOD peaks path for later
    peaks_path = os.path.join(two_nifti, "fod_peaks.nii.gz")

    # Helper function for building paths
    def one_path(subpath):
        return os.path.join(one_raw, subpath)

    # 1. Convert scans to NIFTI or MIF
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

    # Another helper function for building paths
    def five_path(subpath):
        return os.path.join(five_dwi, subpath)

    # 2. Preprocessing steps
    # Combine AP/PA
    run_cmd([
        "mrcat", "-nthreads", str(nthreads),
        five_path("dwi_ap.mif"), five_path("dwi_pa.mif"), five_path("dwi_all.mif"),
        "-axis", "3", "-force"
    ])

    # Denoise
    run_cmd([
        "dwidenoise", "-nthreads", str(nthreads),
        five_path("dwi_all.mif"), five_path("dwi_den.mif"),
        "-noise", five_path("noise.mif"), "-force"
    ])

    # Residual
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        five_path("dwi_all.mif"), five_path("dwi_den.mif"),
        "-subtract", five_path("residual.mif"), "-force"
    ])

    # Gibbs unringing
    run_cmd([
        "mrdegibbs", "-nthreads", str(nthreads),
        five_path("dwi_den.mif"), five_path("dwi_den_unr.mif"),
        "-axes", "0,1", "-force"
    ])

    # Motion/distortion correction with FSL
    run_cmd([
        "dwifslpreproc", "-nthreads", str(nthreads),
        "-rpe_all", "-pe_dir", "AP",
        five_path("dwi_den_unr.mif"), five_path("dwi_den_unr_pre.mif"), "-force"
    ])

    # Bias field correction
    run_cmd([
        "dwibiascorrect", "-nthreads", str(nthreads),
        "ants",
        five_path("dwi_den_unr_pre.mif"), five_path("dwi_den_unr_pre_unbia.mif"),
        "-bias", five_path("bias.mif"), "-force"
    ])

    # Create mask
    run_cmd([
        "dwi2mask", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif"), "-force"
    ])

    # Skull stripping
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif"),
        "-mult", five_path("dwi_den_unr_pre_unbia_skull.mif"), "-force"
    ])

    # Response estimation
    run_cmd([
        "dwi2response", "-nthreads", str(nthreads),
        "dhollander",
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("gm.txt"), five_path("csf.txt"), "-force"
    ])

    # FOD estimation
    run_cmd([
        "dwi2fod", "-nthreads", str(nthreads),
        "msmt_csd",
        "-mask", five_path("mask.mif"),
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("wm.mif"),
        five_path("gm.txt"), five_path("gm.mif"),
        five_path("csf.txt"), five_path("csf.mif"), "-force"
    ])

    # Intensity normalization
    """
    For Boshra: I researched about this command and it works by 
    1) Estimating a polynomial bias field in the log domain (this corrects for intensity inhomogeneities), 
    2) Adjusts the multi-tissue compartment intensities so that their voxel sum converges toward a constant value.
    3) Down weights outliers.
    It's ran individually on each subject and they all end up on a comparable intensity scale, 
    even without a reference subject. 
    """
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        five_path("wm.mif"), five_path("wm_norm.mif"),
        five_path("gm.mif"), five_path("gm_norm.mif"),
        five_path("csf.mif"), five_path("csf_norm.mif"),
        "-mask", five_path("mask.mif"), "-force"
    ])

    # Generate peaks
    run_cmd([
        "sh2peaks",
        five_path("wm_norm.mif"),
        peaks_path,
        "-force"
    ])

    # Run TractSeg for tract segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", subject_dir,
        "--output_type", "tract_segmentation"
    ])

    # Run TractSeg for endings segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", subject_dir,
        "--output_type", "endings_segmentation"
    ])

    # Run TractSeg for Tract Orientation Maps (TOM)
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        ".o", subject_dir,
        "--output_type", "TOM"
    ])

    # Diffusion tensor and ADC/FA computation
    run_cmd([
        "dwi2tensor", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("tensor.mif")
    ])

    run_cmd([
        "tensor2metric", "-nthreads", str(nthreads),
        five_path("tensor.mif"),
        "-adc", five_path("adc.mif")
    ])

    run_cmd([
        "tensor2metric", "-nthreads", str(nthreads),
        five_path("tensor.mif"),
        "-fa", five_path("fa.mif")
    ])


def main():
    # Parser with description for --help
    parser = argparse.ArgumentParser(
        description="Run preproc pipeline for all subjects in a root directory. WILL OVERWRITE FILES"
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
        help="Number of threads to pass to MRtrix commands. Will attempt to use max available threads - 10, if not possible attempts 4."
    )
    args = parser.parse_args()

    # Creates an alphabetically sorted list of absolute paths to directories under given root. Ignores non-directories
    root = os.path.abspath(args.root)
    subject_dirs = sorted([
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    ])

    # Run the pipeline for each subject
    tract_names = pd.read_csv("tract_name.txt", header=None)[0].tolist()
    for subj_dir in subject_dirs:
        #print(f"\n========= Preprocessing subject: {os.path.basename(subj_dir)} =========")
        #preprocess_subject(subj_dir, args.nthreads)
        print(f"\n========= Track generation and resampling subject: {os.path.basename(subj_dir)} =========")
        process_subject(subj_dir,tract_names)


if __name__ == "__main__":
    main()
