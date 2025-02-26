import os
from helpers.helpers import run_cmd


"""
Module with funcitons corresponding to streamline seeding,
tracks, SIFT and TDIs
"""


def streamline_seeding(paths):
    """
    Perform streamline seeding using 5tt2gmwmi.
    """
    os.makedirs(paths["mat_dir"], exist_ok=True)
    fivett_coreg = os.path.join(paths["mat_dir"], "5tt_coreg.mif")
    output_seed = os.path.join(paths["mat_dir"], "gmwmSeed_coreg.mif")
    run_cmd([
        "5tt2gmwmi", fivett_coreg, output_seed, "-fo"
    ])


def generate_tracks_and_sift(paths, nthreads):
    """
    Generate whole-brain tracks and apply SIFT.
    """
    os.makedirs(paths["tck_dir"], exist_ok=True)
    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio.tck")
    fivett_coreg = os.path.join(paths["mat_dir"], "5tt_coreg.mif")
    output_seed = os.path.join(paths["mat_dir"], "gmwmSeed_coreg.mif")
    run_cmd([
        "tckgen", "-act",
        fivett_coreg, "-backtrack",
        "-seed_gmwmi", output_seed,
        "-select", "10000000",
        os.path.join(paths["five_dwi"], "wm_norm.mif"),
        tckgen_output,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # SIFT2 filtering
    """
    FF: Unlike SIFT1, SIFT2 doesn’t prune streamlines but rather calculates a continuous weight
    for each streamline so instead of removing them, SIFT2 scales each one by a factor
    so that the overall tractogram better reflects the fiber densities.
    """
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")
    mu_file = os.path.join(paths["tck_dir"], "mu.txt")
    coeff_file = os.path.join(paths["tck_dir"], "tck_coeffs.txt")
    run_cmd([
        "tcksift2",
        tckgen_output,
        os.path.join(paths["five_dwi"], "wm_norm.mif"),
        sift2_output,
        "-out_mu", mu_file,
        "-out_coeffs", coeff_file,
        "-act", fivett_coreg,
        "-nthreads", str(nthreads),
        "-force"
    ])


def generate_tdis(paths, nthreads):
    """
    Generate Track Density Images :
      - A high-resolution TDI
      - A T1-aligned, scaled TDI using the mu value from SIFT2
    """

    # Define file paths
    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio.tck")
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")
    t1_file = os.path.join(paths["two_nifti"], "t1.mif")
    mu_file = os.path.join(paths["tck_dir"], "mu.txt")
    tdi_output = os.path.join(paths["tck_dir"], "tdi.mif")
    tdi_t1_output = os.path.join(paths["tck_dir"], "tdi_T1.mif")

    # Generate high-resolution TDI.
    run_cmd([
        "tckmap",
        tckgen_output,
        tdi_output,
        "-tck_weights_in", sift2_output,
        "-vox", "0.25",
        "-datatype", "uint16",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Generate T1-aligned and scaled TDI.
    cmd = (
        f"tckmap {tckgen_output} - -template {t1_file} -precise -nthreads {nthreads} | "
        f"mrcalc - $(cat {mu_file}) -mult {tdi_t1_output} -force"
    )
    run_cmd(["bash", "-c", cmd])

