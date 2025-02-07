import os
from helpers import run_cmd


def register_t1_coreg(paths, nthreads):
    """
    Performs the following registration steps:
      1. Extracts and computes b0 from dMRI
      2. Generates 5t segmentation from raw t1
      3. Registers mean b0 to 5tt
      4. Converts the FLIRT transformation to MRtrix format.
      5. Applies the inverse transform to T1 to get the coreg.
    """

    # Mean b0 from dMRI
    dwi_input = os.path.join(paths["five_dwi"], "dwi_den_unr_pre_unbia.mif")
    b0_temp = os.path.join(paths["five_dwi"], "b0.nii.gz")
    mean_b0 = os.path.join(paths["two_nifti"], "mean_b0.nii.gz")
    run_cmd([
        "dwiextract", dwi_input,
        "-bzero", b0_temp, "-force"
    ])

    run_cmd([
        "mrmath", b0_temp, "mean", mean_b0, "-axis", "3", "-force"
    ])

    run_cmd([
        "rm", b0_temp
    ])

    # 5tt T1
    t1_mif = os.path.join(paths["two_nifti"], "t1.mif")
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
    diff2struct_fsl_mat = os.path.join(paths["mat_dir"], "diff2struct_fsl.mat")
    run_cmd([
        "flirt", "-in", mean_b0,
        "-ref", fivett_nocoreg_nii, "-dof", "6",
        "-omat", diff2struct_fsl_mat
    ])

    diff2struct_mrtrix_txt = os.path.join(paths["mat_dir"], "diff2struct_mrtrix.txt")
    run_cmd([
        "transformconvert",
        diff2struct_fsl_mat,
        mean_b0,
        fivett_nocoreg_nii,
        "flirt_import",
        diff2struct_mrtrix_txt,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Apply the inverse transform to raw T1
    t1_coreg = os.path.join(paths["mat_dir"], "t1_coreg.mif")
    run_cmd([
        "mrtransform",
        t1_mif,
        "-linear", diff2struct_mrtrix_txt,
        "-inverse", t1_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])

if __name__ == "__main__":
    """
    Quick code used during troubleshooting 
    """

    import argparse
    from helpers import get_subject_paths

    parser = argparse.ArgumentParser(
        description="Test the register_t1_coreg function in isolation."
    )
    parser.add_argument(
        "--root",
        default="/media/nas/nikita/test_study",
        help="Path to a subject directory with the expected folder structure."
    )
    parser.add_argument(
        "--nthreads",
        type=int,
        default=102,
        help="Number of threads to use."
    )
    args = parser.parse_args()

    root = os.path.abspath(args.root)
    subject_dirs = sorted([
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    ])

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        register_t1_coreg(paths, args.nthreads)
