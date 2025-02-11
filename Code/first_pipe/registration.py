import os
from helpers import run_cmd, get_args


def register_t1_to_dwi(paths, nthreads):
    """
    Generates structural to diffusion transformation matrix
    Registers T1 to diffusion space
    """

    # Mean b0 from dMRI
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
        "mrconvert", "-strides 1,2,3,4",
        mean_b0, mean_b0_newor, "-force"
    ])

    run_cmd([
        "rm", b0_temp
    ])

    # 5tt T1
    t1_mif = os.path.join(paths["two_nifti"], "t1.mif")
    t1_nii = os.path.join(paths["two_nifti"], "t1.nii.gz")
    fivett_nocoreg_mif = os.path.join(paths["two_nifti"], "5tt_nocoreg.mif")
    run_cmd([
        "5ttgen", "fsl", t1_mif,
        fivett_nocoreg_mif, "-force"
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
        "transformconvert",
        struct2diff_fsl,
        t1_nii,
        mean_b0_newor,
        "flirt_import",
        struct2diff_mrtrix,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Apply the inverse transform to raw T1
    t1_coreg = os.path.join(paths["mat_dir"], "t1_coreg.mif")
    run_cmd([
        "mrtransform",
        t1_mif,
        "-linear", struct2diff_mrtrix,
        "-inverse", t1_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])

if __name__ == "__main__":
    """
    Quick code used during troubleshooting 
    """

    from helpers import get_subject_paths

    args = get_args()

    root = os.path.abspath(args.root)
    subject_dirs = sorted([
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    ])

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        register_t1_to_dwi(paths, args.nthreads)
