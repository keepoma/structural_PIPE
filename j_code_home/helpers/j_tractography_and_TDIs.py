import os
from helpers.j_helpers import run_cmd, get_subject_paths


"""
Module with functions corresponding to streamline seeding,
tracks, SIFT and TDIs
"""


def streamline_seeding(paths):
    """
    Perform streamline seeding using 5tt2gmwmi:
    GM/WM Interface Seeding with a Precomputed Mask.
    This differs from dynamic seeding (using -seed_dynamic and -crop_at_gmwmi)
    mainly because the mask provides a more anatomically-informed seeding
    by computing a mask from segmentation data.
    """

    os.makedirs(paths["mat_dir"], exist_ok=True)
    fivett_coreg = os.path.join(paths["dwi_dir"], "5tt_coreg.mif")
    output_seed = os.path.join(paths["dwi_dir"], "gmwmSeed_coreg.mif")
    run_cmd([
        "5tt2gmwmi", fivett_coreg, output_seed, "-fo"
    ])


def generate_tracks_and_sift(paths, nthreads):
    """
    Generate whole-brain tracks and apply SIFT.
    """
    os.makedirs(paths["tck_dir"], exist_ok=True)
    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio_minmax_restricted.tck")
    fivett_coreg = os.path.join(paths["dwi_dir"], "5tt_coreg.mif")
    output_seed = os.path.join(paths["dwi_dir"], "gmwmSeed_coreg.mif")
    wm_norm = os.path.join(paths["dwi_dir"], "wm_norm.mif")
    if not os.path.exists(tckgen_output):
        run_cmd([
            "tckgen",
            "-algorithm", "iFOD2",
            "-act", fivett_coreg, "-backtrack",
            "-seed_gmwmi", output_seed,
            "-select", "10000000",
            "-minlength", "10",
            "-maxlength", "250",
            wm_norm,
            tckgen_output,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print(f"skipping {tckgen_output} gen, already exists")

    # SIFT2 filtering
    """
    FF: Unlike SIFT1, SIFT2 doesn’t prune streamlines but rather calculates a continuous weight
    for each streamline so instead of removing them, SIFT2 scales each one by a factor
    so that the overall tractogram better reflects the fiber densities.
    """
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")
    mu_file = os.path.join(paths["tck_dir"], "mu.txt")
    coeff_file = os.path.join(paths["tck_dir"], "tck_coeffs.txt")
    if not os.path.exists(sift2_output):
        run_cmd([
            "tcksift2",
            tckgen_output,
            os.path.join(paths["dwi_dir"], "wm_norm.mif"),
            sift2_output,
            "-out_mu", mu_file,
            "-out_coeffs", coeff_file,
            "-act", fivett_coreg,
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print(f"skipping {sift2_output} gen, already exists")


if __name__ == '__main__':
    paths = get_subject_paths("/Users/nikitakaruzin/Desktop/Research/Picht/my_brain/me")
    generate_tracks_and_sift(paths, 4)

