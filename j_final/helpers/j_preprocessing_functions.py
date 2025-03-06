import os
import re
from .j_helpers import run_cmd, get_subject_dirs, get_subject_paths


"""
This script contains general purpose importable preprocessing functions
"""


def convert_scans(paths, nthreads):
    """
    Convert anatomical and diffusion scans into standardized NIfTI or MIF formats.
    """
    # Helper lambda for paths from the one_raw directory
    for filename in os.listdir(paths["anat_dir"]):
        if re.search(r'T1w\.nii\.gz', filename):
            t1_nii = os.path.join(paths["anat_dir"], filename)

    t1_mif = os.path.join(paths["anat_dir"], "t1.mif")
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3",
        t1_nii,
        t1_mif,
        "-force"
    ])

    # Convert dMRI AP scan
    for filename in os.listdir(paths["dwi_dir"]):
        if re.search(r'AP\.nii\.gz', filename):
            dwi_ap = os.path.join(paths["dwi_dir"], filename)

    dwi_ap_mif = os.path.join(paths["dwi_dir"], "dwi_ap.mif")
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        dwi_ap, dwi_ap_mif,
        "-force"
    ])

    # Convert dMRI PA scan
    for filename in os.listdir(paths["dwi_dir"]):
        if re.search(r'PA\.nii\.gz', filename):
            dwi_pa = os.path.join(paths["dwi_dir"], filename)

    dwi_pa_mif = os.path.join(paths["dwi_dir"], "dwi_pa.mif")
    run_cmd([
        "mrconvert", "-nthreads", str(nthreads),
        "-strides", "1,2,3,4",
        dwi_pa, dwi_pa_mif,
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

    dwi_path = lambda subpath: os.path.join(paths["dwi_dir"], subpath)

    # Combine AP and PA scans
    run_cmd([
        "mrcat", "-nthreads", str(nthreads),
        dwi_path("dwi_ap.mif"), dwi_path("dwi_pa.mif"), dwi_path("dwi_all.mif"),
        "-axis", "3", "-force"
    ])

    # Denoise
    run_cmd([
        "dwidenoise", "-nthreads", str(nthreads),
        dwi_path("dwi_all.mif"), dwi_path("dwi_den.mif"),
        "-noise", dwi_path("noise.mif"), "-force"
    ])

    # Calculate residual (all - denoised)
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        dwi_path("dwi_all.mif"), dwi_path("dwi_den.mif"),
        "-subtract", dwi_path("residual.mif"), "-force"
    ])

    # Gibbs unringing
    run_cmd([
        "mrdegibbs", "-nthreads", str(nthreads),
        dwi_path("dwi_den.mif"), dwi_path("dwi_den_unr.mif"),
        "-axes", "0,1", "-force"
    ])

    # Motion/distortion correction using FSL
    run_cmd([
        "dwifslpreproc", "-nthreads", str(nthreads),
        "-rpe_all", "-pe_dir", "AP",
        dwi_path("dwi_den_unr.mif"), dwi_path("dwi_den_unr_pre.mif"), "-force"
    ])

    # Bias field correction
    run_cmd([
        "dwibiascorrect", "-nthreads", str(nthreads),
        "ants",
        dwi_path("dwi_den_unr_pre.mif"), dwi_path("dwi_den_unr_pre_unbia.mif"),
        "-bias", dwi_path("bias.mif"), "-force"
    ])

    # Strides correction
    run_cmd([
        "mrconvert", "-strides", "1,2,3,4",
        dwi_path("dwi_den_unr_pre_unbia.mif"),
        dwi_path("dwi_den_unr_pre_unbia_newor.mif"),
        "-force"
    ])

    # Create a mask from the unbiased image
    run_cmd([
        "dwi2mask", "-nthreads", str(nthreads),
        dwi_path("dwi_den_unr_pre_unbia.mif"), dwi_path("mask.mif"), "-force"
    ])

    # Skull stripping (multiplying the image by its mask)
    run_cmd([
        "mrcalc", "-nthreads", str(nthreads),
        dwi_path("dwi_den_unr_pre_unbia.mif"), dwi_path("mask.mif"),
        "-mult", dwi_path("dwi_den_unr_pre_unbia_skull.mif"), "-force"
    ])


def response_function(paths, nthreads):
    """
    Performs response function estimation
    """

    # Helper lambda for constructing file paths
    dwi_path = lambda subpath: os.path.join(paths["dwi_path"], subpath)

    run_cmd([
        "dwi2response", "-nthreads", str(nthreads),
        "dhollander",
        dwi_path("dwi_den_unr_pre_unbia.mif"),
        dwi_path("wm.txt"), dwi_path("gm.txt"), dwi_path("csf.txt"),
        "-force"
    ])


def compute_group_response_functions(root, output_dir, nthreads):
    """
    Compute group-average response functions for each tissue type.
    """

    os.makedirs(output_dir, exist_ok=True)
    subject_dirs = get_subject_dirs(root)

    # List of tissue types to process.
    tissue_types = ["wm", "gm", "csf"]
    response_files = {tissue: [] for tissue in tissue_types}

    # Gather response function files for each subject.
    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            for tissue in tissue_types:
                tissue_file = os.path.join(paths["dwi_dir"], f"{tissue}.txt")
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


def FOD_normalization_peaks(paths, root, nthreads):
    """
    Calculates FOD, performs intensity normalization and generates peaks
    based on individual AND group response functions.
    """

    # Helper lambda for constructing file paths
    dwi_dir = lambda subpath: os.path.join(paths["dwi_dir"], subpath)

    # FOD based on individual RF
    run_cmd([
        "dwi2fod", "-nthreads", str(nthreads),
        "msmt_csd",
        "-mask", dwi_dir("mask.mif"),
        dwi_dir("dwi_den_unr_pre_unbia.mif"),
        dwi_dir("wm.txt"), dwi_dir("wm.mif"),
        dwi_dir("gm.txt"), dwi_dir("gm.mif"),
        dwi_dir("csf.txt"), dwi_dir("csf.mif"),
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
        "-mask", dwi_dir("mask.mif"),
        dwi_dir("dwi_den_unr_pre_unbia.mif"),
        wm_response_file, dwi_dir("wm_group_average_based.mif"),
        gm_response_file, dwi_dir("gm_group_average_based.mif"),
        csf_response_file, dwi_dir("csf_group_average_based.mif"),
        "-force"
    ])

    # Performs global intensity normalization based on individual RF
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        dwi_dir("wm.mif"), dwi_dir("wm_norm.mif"),
        dwi_dir("gm.mif"), dwi_dir("gm_norm.mif"),
        dwi_dir("csf.mif"), dwi_dir("csf_norm.mif"),
        "-mask", dwi_dir("mask.mif"),
        "-force"
    ])

    # Performs global intensity normalization based on group RF
    run_cmd([
        "mtnormalise", "-nthreads", str(nthreads),
        dwi_dir("wm_group_average_based.mif"), dwi_dir("wm_group_average_based_norm.mif"),
        dwi_dir("gm_group_average_based.mif"), dwi_dir("gm_group_average_based_norm.mif"),
        dwi_dir("csf_group_average_based.mif"), dwi_dir("csf_group_average_based_norm.mif"),
        "-mask", dwi_dir("mask.mif"),
        "-force"
    ])

    # Generate peaks on individual and group FODs
    peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")
    peaks_path_group_RF = os.path.join(paths["dwi_dir"], "fod_peaks_group_RF.nii.gz")

    run_cmd([
        "sh2peaks",
        dwi_dir("wm_norm.mif"),
        peaks_path_individual_RF,
        "-force"
    ])

    run_cmd([
        "sh2peaks",
        dwi_dir("wm_group_average_based_norm.mif"),
        peaks_path_group_RF,
        "-force"
    ])


def calculate_tensors_and_dmri_metrics(paths, nthreads):
    """
    Calculate diffusion tensors and derive dMRI metrics (FA, ADC, AD, RD)
    from preprocessed DWI data.
    """

    dwi_image = os.path.join(paths["dwi_dir"], "dwi_den_unr_pre_unbia.mif")
    mask_image = os.path.join(paths["dwi_dir"], "mask.mif")
    tensors_output = os.path.join(paths["dwi_dir"], "tensors.mif")

    # Calculate the diffusion tensor from the preprocessed DWI data
    run_cmd([
        "dwi2tensor",
        dwi_image,
        "-mask", mask_image,
        tensors_output,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Define outputs for the diffusion MRI metrics.
    fa_output_mif = os.path.join(paths["dwi_dir"], "fa.mif")
    fa_output_nii = os.path.join(paths["dwi_dir"], "fa.nii.gz")
    adc_output_mif = os.path.join(paths["dwi_dir"], "adc.mif")
    ad_output_mif = os.path.join(paths["dwi_dir"], "ad.mif")
    rd_output_mif = os.path.join(paths["dwi_dir"], "rd.mif")

    # Fractional Anisotropy map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-fa", fa_output_mif,
        "-nthreads", str(nthreads),
        "-force"
    ])
    run_cmd([
        "mrconvert",
        fa_output_mif,
        fa_output_nii,
        "-force"
    ])

    # Apparent Diffusion Coefficient map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-adc", adc_output_mif,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Calculate the Axial Diffusivity map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-ad", ad_output_mif,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Calculate the Radial Diffusivity map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-rd", rd_output_mif,
        "-nthreads", str(nthreads),
        "-force"
    ])
