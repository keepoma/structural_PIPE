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
        if re.search(r'T1w?\.nii\.gz', filename):
            t1_nii = os.path.join(paths["anat_dir"], filename)
    if not t1_nii:
        raise FileNotFoundError("No T1w file found in " + paths["anat_dir"])

    t1_mif = os.path.join(paths["anat_dir"], "t1.mif")
    if not os.path.exists(t1_mif):
        run_cmd([
            "mrconvert", "-nthreads", str(nthreads),
            "-strides", "1,2,3",
            t1_nii,
            t1_mif,
            "-force"
        ])
    else:
        print(f"{t1_mif} already exists")

    # Convert dMRI AP scan
    dwi_ap_nii = None
    dwi_ap_bvec = None
    dwi_ap_bval = None
    dwi_ap_json = None

    for filename in os.listdir(paths["dwi_dir"]):
        full_path = os.path.join(paths["dwi_dir"], filename)
        if re.search(r'AP.nii\.gz$', filename):
            dwi_ap_nii = full_path
        elif re.search(r'AP.bvec$', filename):
            dwi_ap_bvec = full_path
        elif re.search(r'AP.bval$', filename):
            dwi_ap_bval = full_path
        elif re.search(r'AP.json$', filename):
            dwi_ap_json = full_path
   
    dwi_ap_mif = os.path.join(paths["dwi_dir"], "dwi_ap.mif")
    if not os.path.exists(dwi_ap_mif):

        mrconvert_cmd = [
            "mrconvert", "-nthreads", str(nthreads),
            "-strides", "1,2,3,4",
            dwi_ap_nii, dwi_ap_mif,
            "-force"
        ]
        if dwi_ap_bvec and dwi_ap_bval:
               mrconvert_cmd += ["-fslgrad", dwi_ap_bvec, dwi_ap_bval]
        if dwi_ap_json:
               mrconvert_cmd += ["-json_import", dwi_ap_json]
   
        run_cmd(mrconvert_cmd)
    else:
        print(f"{dwi_ap_mif} already exists")

    # Convert dMRI PA scan
    dwi_pa_nii = None
    dwi_pa_bvec = None
    dwi_pa_bval = None
    dwi_pa_json = None

    for filename in os.listdir(paths["dwi_dir"]):
        full_path = os.path.join(paths["dwi_dir"], filename)
        if re.search(r'PA.nii\.gz$', filename):
            dwi_pa_nii = full_path
        elif re.search(r'PA.bvec$', filename):
            dwi_pa_bvec = full_path
        elif re.search(r'PA.bval$', filename):
            dwi_pa_bval = full_path
        elif re.search(r'PA.json$', filename):
            dwi_pa_json = full_path

    dwi_pa_mif = os.path.join(paths["dwi_dir"], "dwi_pa.mif")
    if not os.path.exists(dwi_pa_mif):

        mrconvert_cmd = [
            "mrconvert", "-nthreads", str(nthreads),
            "-strides", "1,2,3,4",
            dwi_pa_nii,
            dwi_pa_mif,
            "-force"
        ]
        if dwi_pa_bvec and dwi_pa_bval:
            mrconvert_cmd += ["-fslgrad", dwi_pa_bvec, dwi_pa_bval]
        if dwi_pa_json:
            mrconvert_cmd += ["-json_import", dwi_pa_json]

        run_cmd(mrconvert_cmd)
    else:
        print(f"{dwi_pa_mif} already exists")



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

    # Combine AP and PA scans -> output: dwi_all.mif
    if not os.path.exists(dwi_path("dwi_all.mif")):
        run_cmd([
            "mrcat", "-nthreads", str(nthreads),
            dwi_path("dwi_ap.mif"), dwi_path("dwi_pa.mif"), dwi_path("dwi_all.mif"),
            "-axis", "3", "-force"
        ])
    else:
        print("dwi_all.mif already exists. Skipping combining AP/PA scans.")

    # Denoise -> output: dwi_den.mif (and noise.mif)
    if not os.path.exists(dwi_path("dwi_den.mif")):
        run_cmd([
            "dwidenoise", "-nthreads", str(nthreads),
            dwi_path("dwi_all.mif"), dwi_path("dwi_den.mif"),
            "-noise", dwi_path("noise.mif"), "-force"
        ])
    else:
        print("dwi_den.mif already exists. Skipping denoising.")

    # Calculate residual (all - denoised) -> output: residual.mif
    if not os.path.exists(dwi_path("residual.mif")):
        run_cmd([
            "mrcalc", "-nthreads", str(nthreads),
            dwi_path("dwi_all.mif"), dwi_path("dwi_den.mif"),
            "-subtract", dwi_path("residual.mif"), "-force"
        ])
    else:
        print("residual.mif already exists. Skipping residual calculation.")

    # Gibbs unringing -> output: dwi_den_unr.mif
    if not os.path.exists(dwi_path("dwi_den_unr.mif")):
        run_cmd([
            "mrdegibbs", "-nthreads", str(nthreads),
            dwi_path("dwi_den.mif"), dwi_path("dwi_den_unr.mif"),
            "-axes", "0,1", "-force"
        ])
    else:
        print("dwi_den_unr.mif already exists. Skipping Gibbs unringing.")

    # Motion/distortion correction using FSL -> output: dwi_den_unr_pre.mif
    if not os.path.exists(dwi_path("dwi_den_unr_pre.mif")):
        run_cmd([
            "dwifslpreproc", "-nthreads", str(nthreads),
            "-rpe_all", "-pe_dir", "AP",
            dwi_path("dwi_den_unr.mif"), dwi_path("dwi_den_unr_pre.mif"), "-force"
        ])
    else:
        print("dwi_den_unr_pre.mif already exists. Skipping motion/distortion correction.")

    # Bias field correction -> outputs: dwi_den_unr_pre_unbia.mif and bias.mif
    if not (os.path.exists(dwi_path("dwi_den_unr_pre_unbia.mif")) and os.path.exists(dwi_path("bias.mif"))):
        run_cmd([
            "dwibiascorrect", "-nthreads", str(nthreads),
            "ants",
            dwi_path("dwi_den_unr_pre.mif"), dwi_path("dwi_den_unr_pre_unbia.mif"),
            "-bias", dwi_path("bias.mif"), "-force"
        ])
    else:
        print("dwi_den_unr_pre_unbia.mif and bias.mif already exist. Skipping bias field correction.")

    # Strides correction -> output: dwi_den_unr_pre_unbia_newor.mif
    if not os.path.exists(dwi_path("dwi_den_unr_pre_unbia_newor.mif")):
        run_cmd([
            "mrconvert", "-strides", "1,2,3,4",
            dwi_path("dwi_den_unr_pre_unbia.mif"),
            dwi_path("dwi_den_unr_pre_unbia_newor.mif"),
            "-force"
        ])
    else:
        print("dwi_den_unr_pre_unbia_newor.mif already exists. Skipping strides correction.")

    # Create a mask from the unbiased image -> output: mask.mif
    if not os.path.exists(dwi_path("mask.mif")):
        run_cmd([
            "dwi2mask", "-nthreads", str(nthreads),
            dwi_path("dwi_den_unr_pre_unbia.mif"), dwi_path("mask.mif"), "-force"
        ])
    else:
        print("mask.mif already exists. Skipping mask creation.")

    # Skull stripping (multiplying the image by its mask) -> output: dwi_den_unr_pre_unbia_skull.mif
    if not os.path.exists(dwi_path("dwi_den_unr_pre_unbia_skull.mif")):
        run_cmd([
            "mrcalc", "-nthreads", str(nthreads),
            dwi_path("dwi_den_unr_pre_unbia.mif"), dwi_path("mask.mif"),
            "-mult", dwi_path("dwi_den_unr_pre_unbia_skull.mif"), "-force"
        ])
    else:
        print("dwi_den_unr_pre_unbia_skull.mif already exists. Skipping skull stripping.")


def response_function(paths, nthreads):
    """
    Performs response function estimation
    """

    # Helper lambda for constructing file paths
    dwi_path = lambda subpath: os.path.join(paths["dwi_dir"], subpath)
    if not os.path.exists(dwi_path("wm.txt")):
        run_cmd([
            "dwi2response", "-nthreads", str(nthreads),
            "dhollander",
            dwi_path("dwi_den_unr_pre_unbia.mif"),
            dwi_path("wm.txt"), dwi_path("gm.txt"), dwi_path("csf.txt"),
            "-force"
        ])
    else:
        print("wm.txt exists, skipping calc")


def FOD_normalization_peaks(paths, nthreads):
    """
    Calculates FOD, performs intensity normalization and generates peaks
    based on individual AND group response functions.
    """

    # Helper lambda for constructing file paths
    dwi_dir = lambda subpath: os.path.join(paths["dwi_dir"], subpath)

    if not os.path.exists(dwi_dir("fod_peaks_individual_RF.nii.gz")):
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


        # Performs global intensity normalization based on individual RF
        run_cmd([
            "mtnormalise", "-nthreads", str(nthreads),
            dwi_dir("wm.mif"), dwi_dir("wm_norm.mif"),
            dwi_dir("gm.mif"), dwi_dir("gm_norm.mif"),
            dwi_dir("csf.mif"), dwi_dir("csf_norm.mif"),
            "-mask", dwi_dir("mask.mif"),
            "-force"
        ])

        # Generate peaks on individual FODs
        peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")

        run_cmd([
            "sh2peaks",
            dwi_dir("wm_norm.mif"),
            peaks_path_individual_RF,
            "-force"
        ])
    else:
        print("fod_peaks_individual_RF.nii.gz already exists. Skipping mask creation.")


def calculate_tensors_and_dmri_metrics(paths, nthreads):
    """
    Calculate diffusion tensors and derive dMRI metrics (FA, ADC, AD, RD)
    from preprocessed DWI data.
    """

    dwi_image = os.path.join(paths["dwi_dir"], "dwi_den_unr_pre_unbia.mif")
    mask_image = os.path.join(paths["dwi_dir"], "mask.mif")
    tensors_output = os.path.join(paths["dwi_dir"], "tensors.mif")

    # Calculate the diffusion tensor from the preprocessed DWI data
    if not os.path.exists(tensors_output):
        run_cmd([
            "dwi2tensor",
            dwi_image,
            "-mask", mask_image,
            tensors_output,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print(f"{tensors_output} already exists, skipping")

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
