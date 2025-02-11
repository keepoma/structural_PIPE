import os
import pandas as pd
from helpers import run_cmd, get_subject_paths, get_args, prompt_for_folder
from bundle_pipe import process_subject
from registration import register_t1_to_dwi

"""
Main pipeline code. Imports functions from other modules for reusability.
"""


def convert_scans(paths, nthreads, t1_folder, t2_folder, t2_df_folder, dwi_ap_folder, dwi_pa_folder):
    """
    Convert anatomical and diffusion scans into standardized NIfTI or MIF formats.
    """

    # Helper lambda for paths from the one_raw directory
    one_path = lambda subpath: os.path.join(paths["one_raw"], subpath)

    # Convert T1 scan
    t1_nii = os.path.join(paths["two_nifti"], "t1.nii.gz")
    t1_mif = os.path.join(paths["two_nifti"], "t1.mif")
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path(t1_folder), t1_nii,
        "-force"
    ])
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        t1_nii,
        t1_mif,
        "-force"
    ])

    # Convert T2 scan
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path(t2_folder  ), os.path.join(paths["two_nifti"], "t2.nii.gz"),
        "-force"
    ])
    # Convert dark-fluid T2 scan
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        one_path(t2_df_folder),
        os.path.join(paths["two_nifti"], "t2_df.nii.gz"),
        "-force"
    ])
    # Convert dMRI AP scan
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        one_path(dwi_ap_folder), os.path.join(paths["five_dwi"], "dwi_ap.mif"),
        "-force"
    ])
    # Convert dMRI PA scan
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        one_path(dwi_pa_folder), os.path.join(paths["five_dwi"], "dwi_pa.mif"),
        "-force"
    ])


def preprocess_dwi(paths, nthreads):
    """
    Preprocess the dMRI data from the 5_dwi folder. This function includes:
      - Combining AP/PA scans
      - Denoising
      - Calculating residuals
      - Gibbs unringing
      - Motion/distortion correction
      - Bias field correction
      - Mask creation
      - Skull stripping
    """

    # Helper lambda for paths in the 5_dwi folder
    five_path = lambda subpath: os.path.join(paths["five_dwi"], subpath)

    # Combine AP and PA scans
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

    # Calculate residual (all - denoised)
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

    # Motion/distortion correction using FSL
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

    # Strides correction
    run_cmd([
        "mrconvert", "-strides 1,2,3,4",
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("dwi_den_unr_pre_unbia_newor.mif"),
        "-force"
    ])

    # Create a mask from the unbiased image
    run_cmd([
        "dwi2mask", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif"), "-force"
    ])

    # Skull stripping (multiplying the image by its mask)
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"), five_path("mask.mif"),
        "-mult", five_path("dwi_den_unr_pre_unbia_skull.mif"), "-force"
    ])


def fiber_orientation_distribution(paths, nthreads):
    """
    Performs response estimation, FOD estimation, and intensity normalization in one function.
    """

    # Helper lambda for constructing file paths
    five_path = lambda subpath: os.path.join(paths["five_dwi"], subpath)

    run_cmd([
        "dwi2response", "-nthreads", str(nthreads),
        "dhollander",
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("gm.txt"), five_path("csf.txt"),
        "-force"
    ])

    run_cmd([
        "dwi2fod", "-nthreads", str(nthreads),
        "msmt_csd",
        "-mask", five_path("mask.mif"),
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("wm.txt"), five_path("wm.mif"),
        five_path("gm.txt"), five_path("gm.mif"),
        five_path("csf.txt"), five_path("csf.mif"),
        "-force"
    ])

    # Performs global intensity normalization
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        five_path("wm.mif"), five_path("wm_norm.mif"),
        five_path("gm.mif"), five_path("gm_norm.mif"),
        five_path("csf.mif"), five_path("csf_norm.mif"),
        "-mask", five_path("mask.mif"),
        "-force"
    ])

    # Generate peaks
    peaks_path = os.path.join(paths["two_nifti"], "fod_peaks.nii.gz")
    run_cmd([
        "sh2peaks",
        five_path("wm_norm.mif"),
        peaks_path,
        "-force"
    ])


def tractseg_and_tensor(paths, subject_dir, nthreads):
    """
    Performs postprocessing steps for tractography:
    1. Generates FOD peaks from the normalized WM FOD image.
    2. Runs tract segmentation, endings segmentation, and generates tract orientation maps.
    3. Computes the diffusion tensor and derives ADC and FA maps.
    """

    # Helper lambda for paths in the 5_dwi folder
    five_path = lambda subpath: os.path.join(paths["five_dwi"], subpath)
    peaks_path = os.path.join(paths["two_nifti"], "fod_peaks.nii.gz")

    # Tract segmentation
    output_dir = os.path.join(subject_dir, "tractseg_output")
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
        "--output_type", "tract_segmentation"
    ])

    # Endings segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
        "--output_type", "endings_segmentation"
    ])

    #  Tract Orientation Maps
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
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
    # Runs the helper module for cmd arguments
    args = get_args()

    # Creates an alphabetically sorted list of absolute paths to directories under given root. Ignores non-directories
    root = os.path.abspath(args.root)
    subject_dirs = sorted([
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    ])

    # Build the full path to the tract_name.txt file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tract_names_file = os.path.join(script_dir, "tract_name.txt")
    tract_names = pd.read_csv(tract_names_file, header=None)[0].tolist()

    for subj_dir in subject_dirs:
        # Retrieve standard paths
        paths = get_subject_paths(subj_dir)

        # Confirm presence of 1_raw directory
        if not os.path.isdir(paths["one_raw"]):
            print(f"Warning: {paths['one_raw']} does not exist! Exiting.")
            raise SystemExit(1)

        # Creating necessary directories
        os.makedirs(paths["two_nifti"], exist_ok=True)
        os.makedirs(paths["five_dwi"], exist_ok=True)
        os.makedirs(paths["mat_dir"], exist_ok=True)

        # Prompt the user for each folder name (with defaults provided).
        t1_folder = prompt_for_folder("006_T1w_MPR", "T1 scan")
        t2_folder = prompt_for_folder("008_T2w_SPC", "T2 scan")
        t2_df_folder = prompt_for_folder("009_t2_space_dark-fluid_sag_p2_iso_0_8", "dark-fluid T2 scan")
        dwi_ap_folder = prompt_for_folder("016_dMRI_dir98_AP", "dMRI AP scan")
        dwi_pa_folder = prompt_for_folder("019_dMRI_dir98_PA", "dMRI PA scan")

        # Call function to convert scans
        print(f"\n========= Converting Scans for Subject: {os.path.basename(subj_dir)} =========\n")
        convert_scans(paths, args.nthreads, t1_folder, t2_folder, t2_df_folder, dwi_ap_folder, dwi_pa_folder)

        # Call function to preprocess dMRI data
        print(f"\n========= Preprocessing dMRI Data for Subject: {os.path.basename(subj_dir)} =========\n")
        preprocess_dwi(paths, args.nthreads)

        # Call function for fiber orientation distribution
        print(f"\n========= Calculating FOD for Subject: {os.path.basename(subj_dir)} =========\n")
        fiber_orientation_distribution(paths, args.nthreads)

        # Call for tractography postproc
        print(f"\n========= Running tractography for Subject: {os.path.basename(subj_dir)} =========\n")
        tractseg_and_tensor(paths, subj_dir, args.nthreads)

        # Registration of T1 to dMRI using FSL
        print(f"\n========= Registering T1 to dMRI Space for Subject: {os.path.basename(subj_dir)} =========\n")
        register_t1_to_dwi(paths, args.nthreads)

        print(f"\n========= Track generation and resampling Subject: {os.path.basename(subj_dir)} =========\n")
        process_subject(subj_dir, tract_names)
        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")


if __name__ == "__main__":
    main()
