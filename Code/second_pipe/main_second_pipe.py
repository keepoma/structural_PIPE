import os
import Code.preprocess_MRI_data as preproc
import Code.statistical_analysis as sa
from Code.helpers import run_cmd, get_subject_dirs, get_subject_paths, get_args, ask_yes_no, fancy_print
from Code.registration import register_t1_and_5tt_to_dwi


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
        "-force"
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
    atlas_dir = paths["atlas_dir"]
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

    # The labelconvert path will be an issue on non conda mrtrix installations
    parcels_nocoreg = os.path.join(atlas_dir, "hcpmmp1_parcels_nocoreg.mif")
    run_cmd([
        "labelconvert",
        atlas_mif,
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt",
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt",
        parcels_nocoreg
    ])

    # Watchout for errors at this stage, removed "-inverse" flag, implemented struct2diff matrix
    parcels_coreg = os.path.join(atlas_dir, "hcpmmp1_parcels_coreg.mif")
    struct2diff_fsl = os.path.join(paths["mat_dir"], "struct2diff_fsl.mat")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,  # using the parcels generated above
        "-linear", struct2diff_fsl,
        "-datatype", "uint32",
        parcels_coreg
    ])


def connectome_generation(paths):
    """
    Generate the connectome matrix from the filtered tractogram and atlas parcels.
    """

    parcels_coreg = os.path.join(paths["atlas_dir"], "hcpmmp1_parcels_coreg.mif")
    connectome_csv = os.path.join(paths["atlas_dir"], "hcpmmp1.csv")
    assignments_csv = os.path.join(paths["atlas_dir"], "assignments_hcpmmp1.csv")
    run_cmd([
        "tck2connectome",
        "-symmetric",
        "-zero_diagonal",
        "-scale_invnodevol",
        os.path.join(paths["tck_dir"], "sift_1mio.tck"),
        parcels_coreg,
        connectome_csv,
        "-out_assignment", assignments_csv
    ])


def main():
    """
    Main pipeline function: parses arguments, builds subject paths, registers T1 to DWI,
    and processes each subject.
    """

    args = get_args()
    root = os.path.abspath(args.root)
    subject_dirs = get_subject_dirs(root, exclude="group_level_analysis")

    is_preprocessed = ask_yes_no("Is every subject in this folder preprocessed?")
    has_registration = ask_yes_no("Has the registration of T1 and 5tt to dwi been done?")

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        os.makedirs(paths["tck_dir"], exist_ok=True)
        os.makedirs(paths["atlas_dir"], exist_ok=True)
        subject_id = os.path.basename(paths["subject_dir"])

        print(f"\n========= Executing script for Subject: {os.path.basename(subj_dir)} =========")

        if not is_preprocessed:
            fancy_print("Preprocessing", subj_dir)
            fancy_print("Converting Scans", subj_dir)
            preproc.convert_scans(paths, args.nthreads)
            fancy_print("Preprocessing dMRI Data", subj_dir)
            preproc.preprocess_dwi(paths, args.nthreads)
            fancy_print("Calculating Response Function", subj_dir)
            preproc.response_function(paths, args.nthreads)

    group_output_directory = os.path.join(root, "group_analysis")
    print(f"\n========= Calculating Group Response Function =========\n")
    sa.compute_group_response_functions(root, group_output_directory, args.nthreads)

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)

        if not is_preprocessed:
            fancy_print("Performing FOD and normalization", subj_dir)
            preproc.FOD_normalization_peaks(paths, root, args.nthreads)

        if not has_registration:
            fancy_print("Registering T1 to dMRI Space", subj_dir)
            register_t1_and_5tt_to_dwi(paths, args.nthreads)

        fancy_print("Performing streamline seeding", subj_dir)
        streamline_seeding(paths)
        fancy_print("Generating whole-brain tracks and applying SIFT", subj_dir)
        generate_tracks_and_sift(paths, args.nthreads)
        # roi_localization(paths, args.nthreads)
        fancy_print("Generating Freesurfer/HCP-based atlas", subj_dir)
        atlas_generation(paths, args.nthreads, subject_id)
        fancy_print("Generating connectome matrix", subj_dir)
        connectome_generation(paths)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")


if __name__ == "__main__":
    main()
