import os
from .j_helpers import run_cmd


def register_t1_to_dwi(paths, nthreads):
    """
    Generates structural to diffusion transformation matrix
    Registers T1 and 5tt to diffusion space
    """

    # Mean b0 from dMRI
    dwi_input = os.path.join(paths["dwi_dir"], "dwi_den_unr_pre_unbia.mif")
    b0_temp = os.path.join(paths["dwi_dir"], "b0.nii.gz")
    mean_b0 = os.path.join(paths["dwi_dir"], "mean_b0.nii.gz")
    mean_b0_newor = os.path.join(paths["dwi_dir"], "mean_b0_newor.nii.gz")
    t1_nii = os.path.join(paths["anat_dir"], "t1.nii.gz")
    t1_mif = os.path.join(paths["anat_dir"], "t1.mif")


    run_cmd([
        "dwiextract", dwi_input,
        "-bzero", b0_temp, "-force"
    ])

    run_cmd([
        "mrmath", b0_temp, "mean", mean_b0, "-axis", "3", "-force"
    ])

    run_cmd([
        "mrconvert", "-strides", "1,2,3",
        mean_b0, mean_b0_newor, "-force"
    ])

    run_cmd([
        "rm", b0_temp
    ])

    # Matrix with flirt
    struct2diff_fsl = os.path.join(paths["mat_dir"], "struct2diff_fsl.mat")
    os.makedirs(paths["mat_dir"], exist_ok=True)

    run_cmd([
        "flirt", "-in", t1_nii,
        "-ref", mean_b0_newor, "-cost", "normmi",
        "-dof", "6",
        "-omat", struct2diff_fsl
    ])

    struct2diff_mrtrix = os.path.join(paths["mat_dir"], "struct2diff_mrtrix.txt")
    run_cmd([
        "transformconvert", struct2diff_fsl,
        t1_nii, mean_b0_newor,
        "flirt_import", struct2diff_mrtrix,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Apply the transform to raw T1
    t1_coreg = os.path.join(paths["anat_dir"], "t1_coreg.mif")
    run_cmd([
        "mrtransform",
        t1_mif,
        "-linear", struct2diff_mrtrix,
        "-inverse", t1_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])
