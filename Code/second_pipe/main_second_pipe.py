import os
from first_pipe.helpers import run_cmd, get_subject_dirs, get_subject_paths, get_args
from first_pipe.registration import register_t1_and_5tt_to_dwi


def streamline_seeding(paths):
    """
    Perform streamline seeding using 5tt2gmwmi.
    """

    fivett_coreg = os.path.join(paths["mat_dir"], "5tt_coreg.mif")
    output_seed = os.path.join(paths["mat_dir"], "gmwmSeed_coreg.mif")
    run_cmd([
        "5tt2gmwmi", fivett_coreg, output_seed, "-fo"
    ])


def generate_tracks_and_sift(paths, nthreads):
    """
    Generate whole-brain tracks and apply SIFT.
    """

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

    # Create a smaller track file
    smaller_tracks = os.path.join(paths["tck_dir"], "smallerTracks_1000k.tck")
    run_cmd([
        "tckedit", tckgen_output,
        "-number", "1000k",
        smaller_tracks,
        "-nthreads", str(nthreads),
        "-fo"
    ])

    # SIFT filtering on the 10mio tracks
    sift_output = os.path.join(paths["tck_dir"], "sift_1mio.tck")
    run_cmd([
        "tcksift", "-act", fivett_coreg,
        "-term_number", "1000000",
        tckgen_output, os.path.join(paths["five_dwi"], "wm_norm.mif"),
        sift_output,
        "-nthreads", str(nthreads),
        "-fo"
    ])

    # Create a smaller sift file
    smaller_sift = os.path.join(paths["tck_dir"], "smallerSIFT_1000k.tck")
    run_cmd([
        "tckedit", sift_output,
        "-number", "1000k",
        smaller_sift,
        "-nthreads", str(nthreads),
        "-fo"
    ])


def roi_localization(paths, nthreads):
    """
    This funciton is customizable for a specific ROI. Experimental.
    Excluded from main pipeline by default
    """

    # Extract fibers that pass through the ROI
    roi_output = os.path.join(paths["eight_tck"], "FA.tck")
    run_cmd([
        "tckedit",
        "-include", "-26.5,33.95,27.41,3",
        os.path.join(paths["eight_tck"], "sift_1mio.tck"),
        roi_output,
        "-nthreads", str(nthreads),
        "-fo"
    ])


def atlas_generation(paths, nthreads, subject_id):
    """
    Generate the Freesurfer/HCP-based atlas and perform label conversion.
    """

    # Run Freesurfer's recon-all using the subject ID
    t1_nii = os.path.join(paths["two_nifti"], "t1.nii.gz")
    run_cmd([
        "recon-all",
        "-s", subject_id,
        "-i", t1_nii,
        "-all",
        "-threads", str(nthreads)
    ])

    # Map HCP MMP 1.0 atlas from fsaverage onto the subject
    subjects_dir = os.environ.get("SUBJECTS_DIR", "")
    lh_src = os.path.join(subjects_dir, "fsaverage", "label", "lh.hcpmmp1.annot")
    lh_trg = os.path.join(subjects_dir, subject_id, "label", "lh.hcpmmp1.annot")
    rh_src = os.path.join(subjects_dir, "fsaverage", "label", "rh.hcpmmp1.annot")
    rh_trg = os.path.join(subjects_dir, subject_id, "label", "rh.hcpmmp1.annot")
    run_cmd([
        "mri_surf2surf",
        "--srcsubject", "fsaverage",
        "--trgsubject", subject_id,
        "--hemi", "lh",
        "--sval-annot", lh_src,
        "--tval", lh_trg
    ])
    run_cmd([
        "mri_surf2surf",
        "--srcsubject", "fsaverage",
        "--trgsubject", subject_id,
        "--hemi", "rh",
        "--sval-annot", rh_src,
        "--tval", rh_trg
    ])

    # Create atlas outputs in the designated raw atlas folder
    atlas_dir = paths["raw_atlas"]
    os.makedirs(atlas_dir, exist_ok=True)
    atlas_mgz = os.path.join(atlas_dir, "hcpmmp1.mgz")
    run_cmd([
        "mri_aparc2aseg",
        "--old-ribbon",
        "-s", subject_id,
        "--annot", "hcpmmp1",
        "-o", atlas_mgz
    ])
    atlas_mif = os.path.join(atlas_dir, "hcpmmp1.mif")
    run_cmd([
        "mrconvert",
        "-datatype", "uint32",
        atlas_mgz,
        atlas_mif
    ])
    parcels_nocoreg = os.path.join(atlas_dir, "hcpmmp1_parcels_nocoreg.mif")
    run_cmd([
        "labelconvert",
        atlas_mif,
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt",
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt",
        parcels_nocoreg
    ])
    parcels_coreg = os.path.join(atlas_dir, "hcpmmp1_parcels_coreg.mif")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,  # using the parcels generated above
        "-linear", transform_mat,
        "-inverse",
        "-datatype", "uint32",
        parcels_coreg
    ])


def connectome_generation(paths, nthreads):
    """
    Generate the connectome matrix from the filtered tractogram and atlas parcels.
    """
    parcels_coreg = os.path.join(paths["raw_atlas"], "hcpmmp1_parcels_coreg.mif")
    connectome_csv = os.path.join(paths["raw_atlas"], "hcpmmp1.csv")
    assignments_csv = os.path.join(paths["raw_atlas"], "assignments_hcpmmp1.csv")
    run_cmd([
        "tck2connectome",
        "-symmetric",
        "-zero_diagonal",
        "-scale_invnodevol",
        os.path.join(paths["eight_tck"], "sift_1mio.tck"),
        parcels_coreg,
        connectome_csv,
        "-out_assignment", assignments_csv
    ])


def process_subject(paths, nthreads):
    """
    Process an individual subject through seeding, tracking, ROI localization,
    atlas generation, and connectome matrix generation.
    """
    # Do not change the working directory; instead, work with absolute paths.
    subject_id = os.path.basename(paths["subject_dir"])

    streamline_seeding(paths, nthreads)
    generate_tracks_and_sift(paths, nthreads)
    #roi_localization(paths, nthreads)
    atlas_generation(paths, nthreads, subject_id)
    connectome_generation(paths, nthreads)


def main():
    """
    Main pipeline function: parses arguments, builds subject paths, registers T1 to DWI,
    and processes each subject.
    """
    args = get_args()
    root = os.path.abspath(args.root)
    subject_dirs = get_subject_dirs(root, exclude="group_level_analysis")

    for subj_dir in subject_dirs:
        print(f"\n========= Processing Subject: {os.path.basename(subj_dir)} =========\n")

        # Retrieve standard paths and update with additional folder keys
        paths = get_subject_paths(subj_dir)
        paths.update({
            "subject_dir": subj_dir,
            "six_mif": os.path.join(subj_dir, "IN", "6_mif"),
            "eight_tck": os.path.join(subj_dir, "IN", "8_tck"),
            "five_dwi": os.path.join(subj_dir, "IN", "5_dwi"),
            "seven_mat": os.path.join(subj_dir, "IN", "7_mat"),
            "two_nifti": os.path.join(subj_dir, "IN", "2_nifti"),
            "raw_atlas": os.path.join(subj_dir, "raw", "9_atlas")
        })

        # Create directories if they do not exist
        for key in ["six_mif", "eight_tck", "five_dwi", "seven_mat", "two_nifti", "raw_atlas"]:
            os.makedirs(paths[key], exist_ok=True)

        # Register T1 to dMRI using the provided registration function.
        register_t1_and_5tt_to_dwi(paths, args.nthreads)

        # Process the subject through the pipeline.
        process_subject(paths, args.nthreads)


if __name__ == "__main__":
    main()
