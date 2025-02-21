import os
from Code.helpers import run_cmd, get_subject_paths, get_subject_dirs, get_args, prompt_for_folder


"""
This script contains general purpose importable preprocessing functions
"""


def convert_scans(paths, nthreads):
    """
    Convert anatomical and diffusion scans into standardized NIfTI or MIF formats.
    """

    # Pending to save these to config file or something in case of script re-run
    print("Please provide the following folder names: ")
    t1_folder = prompt_for_folder("006_T1w_MPR", "T1 scan")
    t2_folder = prompt_for_folder("008_T2w_SPC", "T2 scan")
    t2_df_folder = prompt_for_folder("009_t2_space_dark-fluid_sag_p2_iso_0_8", "T2 FLAIR")
    dwi_ap_folder = prompt_for_folder("016_dMRI_dir98_AP", "dMRI AP scan")
    dwi_pa_folder = prompt_for_folder("019_dMRI_dir98_PA", "dMRI PA scan")

    # Helper lambda for paths from the one_raw directory
    one_path = lambda subpath: os.path.join(paths["one_raw"], subpath)

    # Convert T1 scan
    os.makedirs(paths["two_nifti"], exist_ok=True)
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
    os.makedirs(paths["five_dwi"], exist_ok=True)
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
        "mrconvert", "-strides", "1,2,3,4",
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


def response_function(paths, nthreads):
    """
    Performs response function estimation
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


def FOD_normalization_peaks(paths, root, nthreads):
    """
    Calculates FOD, performs intensity normalization and generates peaks
    based on individual AND group response functions.
    """

    # Helper lambda for constructing file paths
    five_path = lambda subpath: os.path.join(paths["five_dwi"], subpath)

    # FOD based on individual RF
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

    # FOD based on group RF
    group_output_directory = os.path.join(root, "group_analysis")
    wm_response_file = os.path.join(group_output_directory, "group_average_response_wm.txt")
    gm_response_file = os.path.join(group_output_directory, "group_average_response_gm.txt")
    csf_response_file = os.path.join(group_output_directory, "group_average_response_csf.txt")
    run_cmd([
        "dwi2fod", "-nthreads", str(nthreads),
        "msmt_csd",
        "-mask", five_path("mask.mif"),
        five_path("dwi_den_unr_pre_unbia.mif"),
        wm_response_file, five_path("wm_group_average_based.mif"),
        gm_response_file, five_path("gm_group_average_based.mif"),
        csf_response_file, five_path("csf_group_average_based.mif"),
        "-force"
    ])

    # Performs global intensity normalization based on individual RF
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        five_path("wm.mif"), five_path("wm_norm.mif"),
        five_path("gm.mif"), five_path("gm_norm.mif"),
        five_path("csf.mif"), five_path("csf_norm.mif"),
        "-mask", five_path("mask.mif"),
        "-force"
    ])

    # Performs global intensity normalization based on group RF
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        five_path("wm_group_average_based.mif"), five_path("wm_group_average_based_norm.mif"),
        five_path("gm_group_average_based.mif"), five_path("gm_group_average_based_norm.mif"),
        five_path("csf_group_average_based.mif"), five_path("csf_group_average_based_norm.mif"),
        "-mask", five_path("mask.mif"),
        "-force"
    ])

    # Generate peaks
    peaks_path_individual_RF = os.path.join(paths["two_nifti"], "fod_peaks_individual_RF.nii.gz")
    peaks_path_group_RF = os.path.join(paths["two_nifti"], "fod_peaks_group_RF.nii.gz")

    run_cmd([
        "sh2peaks",
        five_path("wm_norm.mif"),
        peaks_path_individual_RF,
        "-force"
    ])

    run_cmd([
        "sh2peaks",
        five_path("wm_group_average_based_norm.mif"),
        peaks_path_group_RF,
        "-force"
    ])
