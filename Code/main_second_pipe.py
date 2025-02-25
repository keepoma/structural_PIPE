import os
import preprocess_MRI_data as preproc
import statistical_analysis as sa
from helpers.helpers import (run_cmd, get_subject_dirs,
                             get_subject_paths,get_args,
                             ask_yes_no, fancy_print)
from registration import register_t1_and_5tt_to_dwi

"""
Main pipeline for connectome construction
Example run: python3 main_second_pipe.py --root /home/nikita/Nikita_MRI

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
    # Replaced this code with more current tcksift2
    # FF: SIFT1 prunes streamlines. After SIFT1, all remaining streamlines have the same weight (1)
    """
    sift2_output = os.path.join(paths["tck_dir"], "sift_1mio.tck")
    run_cmd([
        "tcksift", "-act", fivett_coreg,
        "-term_number", "1000000",
        tckgen_output, os.path.join(paths["five_dwi"], "wm_norm.mif"),
        sift2_output,
        "-nthreads", str(nthreads),
        "-fo"
    ])
    """

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


def roi_localization(paths, nthreads):
    """
    This function is customizable for a specific ROI.
    Experimental, excluded from main pipeline by default
    """

    # Extract fibers that pass through the custom ROI
    roi_output = os.path.join(paths["eight_tck"], "FA.tck")
    run_cmd([
        "tckedit",
        "-include", "-26.5,33.95,27.41,3",
        os.path.join(paths["eight_tck"], "sift_1mio.tck"),
        roi_output,
        "-nthreads", str(nthreads),
        "-fo"
    ])


def freesurfer_atlas_generation(paths, nthreads, subject_id):
    """
    Generate the Freesurfer/HCP-based atlas and perform label conversion.
    """

    os.makedirs(paths["atlas_dir"], exist_ok=True)
    # Run Freesurfer's recon-all using the subject ID
    # use command rm -rf $SUBJECTS_DIR/subject_id if rerunning
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

    # Create atlas outputs in the designated atlas folder
    atlas_dir = paths["atlas_dir"]
    atlas_mgz = os.path.join(atlas_dir, "hcpmmp1.mgz")
    atlas_mif = os.path.join(atlas_dir, "hcpmmp1.mif")


    run_cmd([
        "mri_aparc2aseg",
        "--old-ribbon",
        "--s", subject_id,
        "--annot", "hcpmmp1",
        "--o", atlas_mgz
    ])

    run_cmd([
        "mrconvert",
        "-datatype", "uint32",
        atlas_mgz,
        atlas_mif,
        "-force"
    ])

    # The labelconvert path will be an issue on non-conda mrtrix installations
    parcels_nocoreg = os.path.join(atlas_dir, "hcpmmp1_parcels_nocoreg.mif")
    run_cmd([
        "labelconvert",
        atlas_mif,
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt",
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt",
        parcels_nocoreg,
        "-force"
    ])

    # Watchout for errors at this stage, removed "-inverse" flag, implemented struct2diff matrix
    parcels_coreg = os.path.join(atlas_dir, "hcpmmp1_parcels_coreg.mif")
    struct2diff_fsl = os.path.join(paths["mat_dir"], "struct2diff_fsl.mat")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,  # using the parcels generated above
        "-linear", struct2diff_fsl,
        "-datatype", "uint32",
        parcels_coreg,
        "-force"
    ])


def nextbrain_atlas_generation(paths, nthreads, subject_id):
    """
    Generate the NextBrain atlas and prepare label conversion outputs.

    Assumed to produce the following key files in paths atlas_dir:
      - seg.left.mgz and seg.right.mgz: NextBrain segmentations for the left and right hemispheres.
      - vols.left.csv and vols.right.csv: ROI volumes from the soft (posterior) segmentations.
      - seg.left.rgb.mgz: a color-coded version of the soft segmentation.
      - synthseg.mgz and SynthSeg_volumes.csv: outputs from the initial SynthSeg run.
      - atlas_reg.nii.gz: registration to MNI space.
      - lut.txt: a lookup table mapping label indices to anatomical names.
      - atlas_nonlinear_reg.[left/right].nii.gz: registered cartoon images (for debugging).

    For connectome generation, this function further merges the left/right segmentations,
    converts the combined volume to MRtrix (.mif) format, performs label conversion using the
    provided LUT, and applies a linear transform (e.g. from structural to diffusion space).
    """

    os.makedirs(paths["atlas_dir"], exist_ok=True)
    t1_nii = os.path.join(paths["two_nifti"], "t1.nii.gz")


    """
    # Run NextBrain segmentation
    # This command is supposed to work with FireANTs or SynthMorph
    # SynthMorph example: 
    # mri_histo_atlas_segment_fast $SUBJECTS_DIR/bert/mri/orig.mgz /path/to/output/bert_histo_atlas_segmentation/ 1 8
    """
    run_cmd([
        "mri_histo_atlas_segment_fast",
        t1_nii, os.path.join(paths["atlas_dir"], "t1.nii.gz"),
        "1", str(nthreads)
    ])

    # Merge left/right hemisphere segmentations
    left_seg = os.path.join(paths["atlas_dir"], "seg.left.mgz")
    right_seg = os.path.join(paths["atlas_dir"], "seg.right.mgz")
    combined_seg = os.path.join(paths["atlas_dir"], "seg.mgz")
    run_cmd([
        "mri_concat",
        "--o", combined_seg,
        left_seg,
        right_seg
    ])

    # Convert combined segmentation to .mif
    combined_seg_mif = os.path.join(paths["atlas_dir"], "nextbrain_seg.mif")
    run_cmd([
        "mrconvert",
        "-datatype", "uint32",
        combined_seg,
        combined_seg_mif,
        "-threads", str(nthreads),
        "-force"
    ])

    # Perform label conversion using the NextBrain lookup table
    # Here  assume that NextBrain outputs a file 'lut.txt' in the atlas directory.
    # If NextBrain does not provide separate original/ordered LUTs, we use the same file for both.
    lut_file = os.path.join(paths["atlas_dir"], "lut.txt")
    parcels_nocoreg = os.path.join(paths["atlas_dir"], "nextbrain_parcels_nocoreg.mif")
    run_cmd([
        "labelconvert",
        combined_seg_mif,
        lut_file,  # source LUT
        lut_file,  # target LUT (assumed to be the same)
        parcels_nocoreg,
        "-threads", str(nthreads),
        "-force"
    ])

    # Coregister the parcels to diffusion space
    struct2diff = os.path.join(paths["mat_dir"], "struct2diff.mat")
    parcels_coreg = os.path.join(paths["atlas_dir"], "nextbrain_parcels_coreg.mif")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,
        "-linear", struct2diff,
        "-datatype", "uint32",
        parcels_coreg,
        "-threads", str(nthreads),
        "-force"
    ])


def connectome_generation(paths, nthreads):
    """
    Generate the connectome matrix from the filtered tractogram and atlas parcels.
    Nodes and edges
    """

    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio.tck")
    parcels_nocoreg = os.path.join(paths["atlas_dir"], "hcpmmp1_parcels_nocoreg.mif")
    connectome_csv = os.path.join(paths["atlas_dir"], "hcpmmp1.csv")
    assignments_csv = os.path.join(paths["atlas_dir"], "assignments_hcpmmp1.csv")
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")

    """
    FF: As far as I understand, we dont use -scale_invnodevol nor -scale_invlenght,
    because f.e. -scale_invnodevol tells tck2connectome to scale each connection by 
    the inverse of the parcel volume, but we don't need that because we have the
    external streamline weights from SIFT2 now
    Are these combinable? 
    """
    # It would be interesting to consider scaling by parcel volume or connection length with
    # -scale_invnodevol or -scale_invlength, respectively
    run_cmd([
        "tck2connectome",
        "-tck_weights_in", sift2_output,
        tckgen_output,
        parcels_nocoreg,
        connectome_csv,
        "-out_assignment", assignments_csv,
        "-symmetric", "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])
    #mrview hcpmmp1_parcels_coreg.mif -connectome.init hcpmmp1_parcels_coreg.mif -connectome.load hcpmmp1.csv

    # Representing nodes as cortical meshes
    mesh_obj = os.path.join(paths["atlas_dir"], "hcpmmp1_mesh.obj")
    run_cmd([
        "label2mesh", parcels_nocoreg, mesh_obj,
        "-force"
    ])

    # Constructing a representation of true edge routes with exemplars
    # I think Boshra uses the nocoreg parcels, make sure to check this step. Changed back to boshras
    exemplars = os.path.join(paths["tck_dir"], "exemplars.tck")
    run_cmd([
        "connectome2tck", tckgen_output, assignments_csv,
        exemplars, "-tck_weights_in", sift2_output,
        "-exemplars", parcels_nocoreg,
        "-files", "single",
        "-nthreads", str(nthreads),
        "-force"
    ])


def calculate_tensors_and_dmri_metrics(paths, nthreads):
    """
    Calculate diffusion tensors and derive dMRI metrics (FA, ADC, AD, RD)
    from preprocessed DWI data.
    """

    dwi_image = os.path.join(paths["five_dwi"], "dwi_den_unr_pre_unbia.mif")
    mask_image = os.path.join(paths["five_dwi"], "mask.mif")
    tensors_output = os.path.join(paths["five_dwi"], "tensors.mif")

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
    fa_output = os.path.join(paths["five_dwi"], "fa.mif")
    adc_output = os.path.join(paths["five_dwi"], "adc.mif")
    ad_output = os.path.join(paths["five_dwi"], "ad.mif")
    rd_output = os.path.join(paths["five_dwi"], "rd.mif")

    # Fractional Anisotropy map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-fa", fa_output,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Apparent Diffusion Coefficient map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-adc", adc_output,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Calculate the Axial Diffusivity map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-ad", ad_output,
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Calculate the Radial Diffusivity map
    run_cmd([
        "tensor2metric",
        tensors_output,
        "-mask", mask_image,
        "-rd", rd_output,
        "-nthreads", str(nthreads),
        "-force"
    ])


def generate_weighted_connectome_matrices(paths, nthreads):
    """
    Generate dMRI metrics along streamlines and construct connectome matrices weighted by these metrics.

    This function performs the following steps:
      1. For each diffusion metric (FA, ADC, AD, RD), sample the underlying image along the streamlines
         to compute a mean value per streamline
      2. Generate connectome matrices weighted by each metric by scaling the streamline contributions
      3. Calculate the average streamline length between node pairs
    """

    tck_file = os.path.join(paths["tck_dir"], "tracks_10mio.tck")
    parcels_file = os.path.join(paths["atlas_dir"], "hcpmmp1_parcels_nocoreg.mif")
    weights_file = os.path.join(paths["tck_dir"], "sift2weights.csv")

    # Define output file paths for the streamline metric sampling.
    os.makedirs(paths["connectome_dir"], exist_ok=True)
    out_fa = os.path.join(paths["connectome_dir"], "tck_meanFA.csv")
    out_adc = os.path.join(paths["connectome_dir"], "tck_meanADC.csv")
    out_ad = os.path.join(paths["connectome_dir"], "tck_meanAD.csv")
    out_rd = os.path.join(paths["connectome_dir"], "tck_meanRD.csv")

    # Sample dMRI metrics along the streamlines (tcksample with "-stat_tck mean").
    run_cmd([
        "tcksample",
        tck_file,
        os.path.join(paths["five_dwi"], "fa.mif"),
        out_fa,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tck_file,
        os.path.join(paths["five_dwi"], "adc.mif"),
        out_adc,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tck_file,
        os.path.join(paths["five_dwi"], "ad.mif"),
        out_ad,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tck_file,
        os.path.join(paths["five_dwi"], "rd.mif"),
        out_rd,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Generate connectome matrices weighted by each metric.
    # For FA-weighted connectome:
    run_cmd([
        "tck2connectome",
        tck_file,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_fa.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_fa,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For ADC-weighted connectome:
    run_cmd([
        "tck2connectome",
        tck_file,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_adc.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_adc,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For AD-weighted connectome:
    run_cmd([
        "tck2connectome",
        tck_file,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_ad.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_ad,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For RD-weighted connectome:
    run_cmd([
        "tck2connectome",
        tck_file,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_rd.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_rd,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Generate a connectome matrix of mean streamline lengths
    run_cmd([
        "tck2connectome",
        tck_file,
        parcels_file,
        os.path.join(paths["connectome_dir"], "meanlength.csv"),
        "-tck_weights_in", weights_file,
        "-scale_length",
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])


def main():
    """
    Main pipeline function: parses arguments, builds subject paths, registers T1 to DWI,
    and processes each subject.
    """

    args = get_args()
    root = os.path.abspath(args.root)
    subject_dirs = get_subject_dirs(root, "group_analysis")

    is_preprocessed = ask_yes_no("Is every subject in this folder preprocessed?")
    has_registration = ask_yes_no("Has the registration of T1 and 5tt to dwi been done?")

    for subj_dir in subject_dirs:
        # Retrieve standard paths
        paths = get_subject_paths(subj_dir)
        subject_id = os.path.basename(subj_dir)
        fancy_print("Executing script for Subject:", subj_dir)

        if not is_preprocessed:
            fancy_print("Preprocessing", subj_dir)
            fancy_print("Converting Scans", subj_dir)
            preproc.convert_scans(paths, args.nthreads)
            fancy_print("Preprocessing dMRI Data", subj_dir)
            preproc.preprocess_dwi(paths, args.nthreads)
            fancy_print("Calculating Response Function", subj_dir)
            preproc.response_function(paths, args.nthreads)

    if len(subject_dirs) > 1:
        group_output_directory = os.path.join(root, "group_analysis")
        print(f"\n========= Calculating Group Response Function =========\n")
        sa.compute_group_response_functions(root, group_output_directory, args.nthreads)

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)

        if not is_preprocessed:
            fancy_print("Performing FOD and normalization", subj_dir)
            preproc.FOD_normalization_peaks(paths, root, args.nthreads)

        if not has_registration:
            fancy_print("Registering T1 and 5tt to dMRI Space", subj_dir)
            register_t1_and_5tt_to_dwi(paths, args.nthreads)

        """
        fancy_print("Performing streamline seeding", subj_dir)
        streamline_seeding(paths)
        fancy_print("Generating whole-brain tracks and applying SIFT", subj_dir)
        generate_tracks_and_sift(paths, args.nthreads)
        fancy_print("Generating TDIs and aligning T1", subj_dir)
        generate_tdis(paths, args.nthreads)
        # roi_localization(paths, args.nthreads)
        fancy_print("Generating Freesurfer/HCP-based atlas", subj_dir)
        freesurfer_atlas_generation(paths, args.nthreads, subject_id)
        fancy_print("Generating connectome matrix", subj_dir)
        connectome_generation(paths, args.nthreads)
        fancy_print("Calculating Tensor and related metrics", subj_dir)
        calculate_tensors_and_dmri_metrics(paths, args.nthreads)
        """
        fancy_print("Generating connectome weighting by metrics", subj_dir)
        generate_weighted_connectome_matrices(paths, args.nthreads)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")


if __name__ == "__main__":
    main()
