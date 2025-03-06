import os
from .j_helpers import run_cmd


def register_t1_and_5tt_to_dwi(paths, nthreads, do_hsvs=False):
    """
    Generates structural to diffusion transformation matrix
    Registers T1 and 5tt to diffusion space
    """

    # Mean b0 from dMRI
    os.makedirs(paths["mat_dir"], exist_ok=True)
    dwi_input = os.path.join(paths["five_dwi"], "dwi_den_unr_pre_unbia.mif")
    b0_temp = os.path.join(paths["five_dwi"], "b0.nii.gz")
    mean_b0 = os.path.join(paths["two_nifti"], "mean_b0.nii.gz")
    mean_b0_newor = os.path.join(paths["two_nifti"], "mean_b0_newor.nii.gz")

    run_cmd([
        "dwiextract", dwi_input,
        "-bzero", b0_temp, "-force"
    ])

    run_cmd([
        "mrmath", b0_temp, "mean", mean_b0, "-axis", "3", "-force"
    ])

    run_cmd([
        "mrconvert", "-strides", "1,2,3,4",
        mean_b0, mean_b0_newor, "-force"
    ])

    run_cmd([
        "rm", b0_temp
    ])

    # 5tt T1 fsl
    t1_mif = os.path.join(paths["two_nifti"], "t1.mif")
    t1_nii = os.path.join(paths["two_nifti"], "t1.nii.gz")
    fivett_nocoreg_mif = os.path.join(paths["two_nifti"], "5tt_nocoreg.mif")
    run_cmd([
        "5ttgen", "fsl", t1_mif,
        fivett_nocoreg_mif, "-force"
    ])

    # 5tt T1 hsvs
    if do_hsvs:
        subjects_dir = os.environ.get("SUBJECTS_DIR", "")
        fivett_nocoreg_hsvs_mif = os.path.join(paths["two_nifti"], "5tt_nocoreg_hsvs.mif")
        run_cmd([
            "5ttgen", "hsvs", subjects_dir,
            fivett_nocoreg_hsvs_mif,
            "-template", t1_mif,
            "-nocrop"
        ])

    fivett_nocoreg_nii = os.path.join(paths["two_nifti"], "5tt_nocoreg.nii.gz")
    run_cmd([
        "mrconvert", fivett_nocoreg_mif, fivett_nocoreg_nii,
        "-force"
    ])

    # Matrix with flirt
    struct2diff_fsl = os.path.join(paths["mat_dir"], "struct2diff_fsl.mat")
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
    t1_coreg = os.path.join(paths["mat_dir"], "t1_coreg.mif")
    run_cmd([
        "mrtransform",
        t1_mif,
        "-linear", struct2diff_mrtrix,
        "-inverse", t1_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Apply the transform to 5tt
    fivett_coreg = os.path.join(paths["mat_dir"], "5tt_coreg.mif")
    run_cmd([
        "mrtransform", fivett_nocoreg_mif,
        "-linear", struct2diff_mrtrix,
        "-inverse", fivett_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])