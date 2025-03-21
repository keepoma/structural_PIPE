import os
from helpers.helpers import run_cmd, get_subject_paths, get_subject_dirs


"""
Functions for atlas generation.
Currently based on:
    -NextBrain
    -Freesurfer
"""


def nextbrain_atlas_generation(paths, nthreads):
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
    subjects_dir = os.environ.get("SUBJECTS_DIR", "")
    orig_mgz = os.path.join(subjects_dir, "me", "mri", "orig.mgz")
    orig_mgz = "/home/nikita/freesurfer/subjects/me/mri/orig.mgz"


    """
    # Run NextBrain segmentation
    # This command is supposed to work with FireANTs or SynthMorph
    # SynthMorph example: 
    # mri_histo_atlas_segment_fast $SUBJECTS_DIR/bert/mri/orig.mgz /path/to/output/bert_histo_atlas_segmentation/ 1 8
    """
    """
    run_cmd([
        "mri_histo_atlas_segment",
        orig_mgz, os.path.join(paths["atlas_dir"], "nextbrain_segmentation"),
        "simplified",
        "1", str(nthreads)
    ])
    """
    
    """
    mri_convert \
  -rl /home/nikita/freesurfer/subjects/me/mri/orig.mgz \
  /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_left.mgz \
  /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_left_in_orig.mgz
  
    mri_convert \
  -rl /home/nikita/freesurfer/subjects/me/mri/orig.mgz \
  /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_right.mgz \
  /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_right_in_orig.mgz  

    """

    # Merge left/right hemisphere segmentations
    left_seg = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "seg_left.mgz")
    right_seg_res = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "seg_right_resampled.mgz")
    combined_seg = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "seg.mgz")

    """
    mri_concat \
  --i /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_left_in_orig.mgz \
  --i /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg_right_in_orig.mgz \
  --o /home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/seg.mgz \
  --max
    """
    """
    run_cmd([
        "mri_concat",
        "--o", combined_seg,
        left_seg,
        right_seg_res
    ])   
    """

    # Convert combined segmentation to .mif
    combined_seg_mif = os.path.join(paths["atlas_dir"], "nextbrain_seg.mif")
    run_cmd([
        "mrconvert",
        "-datatype", "uint32",
        combined_seg,
        combined_seg_mif,
        "-force"
    ])
    
    """
    mrmath /home/nikita/Nikita_MRI/me/atlas/nextbrain_seg.mif -axis 3 max /home/nikita/Nikita_MRI/me/atlas/nextbrain_seg_3D.mif
    """

    # Perform label conversion using the NextBrain lookup table
    # Here  assume that NextBrain outputs a file 'lut.txt' in the atlas directory.
    # If NextBrain does not provide separate original/ordered LUTs, we use the same file for both.
    lut_file = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "lookup_table.txt")
    parcels_nocoreg = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "nextbrain_parcels_nocoreg.mif")
    combined_seg_mif_3d = os.path.join(paths["atlas_dir"], "nextbrain_seg_3D.mif")
    run_cmd([
        "labelconvert",
        combined_seg_mif,
        lut_file,  # source LUT
        lut_file,  # target LUT (assumed to be the same)
        parcels_nocoreg,
        "-force"
    ])

    # Coregister the parcels to diffusion space
    struct2diff_txt = os.path.join(paths["mat_dir"], "struct2diff_mrtrix.txt")
    parcels_coreg = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "nextbrain_parcels_coreg.mif")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,
        "-linear", struct2diff_txt,
        "-datatype", "uint32",
        parcels_coreg,
        "-force"
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
        "-strides", "1,2,3",
        atlas_mif,
        "-force"
    ])


    # The labelconvert path will be an issue on non-conda mrtrix installations
    parcels_nocoreg = os.path.join(atlas_dir, "hcpmmp1_parcels_nocoreg.mif")
    run_cmd([
        "labelconvert",
        atlas_mif,
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt",
        #"/Users/nikitakaruzin/MRI/projects/BATMAN/Supplementary_Files/hcpmmp1_original.txt",
        "/home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt",
        #"/Users/nikitakaruzin/MRI/projects/BATMAN/Supplementary_Files/hcpmmp1_ordered.txt",
        parcels_nocoreg,
        "-force"
    ])   



    # Watchout for errors at this stage, removed "-inverse" flag, implemented struct2diff matrix
    parcels_coreg = os.path.join(atlas_dir, "hcpmmp1_parcels_coreg.mif")
    struct2diff_fsl = os.path.join(paths["mat_dir"], "struct2diff_fsl.mat")
    # Experimental txt file
    struct2diff_txt = os.path.join(paths["mat_dir"], "struct2diff_mrtrix.txt")
    run_cmd([
        "mrtransform",
        parcels_nocoreg,  # using the parcels generated above
        "-linear", struct2diff_txt,
        "-datatype", "uint32",
        parcels_coreg,
        "-force"
    ])


"""
if __name__ == '__main__':
    paths = get_subject_paths("/home/nikita/Nikita_MRI/me")
    nextbrain_atlas_generation(paths, 50)
    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio.tck")
    parcels_coreg = "/home/nikita/Nikita_MRI/me/atlas/nextbrain_segmentation/nextbrain_parcels_coreg.mif"
    connectome_csv = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "nextbrain_connectome.csv")
    assignments_csv = os.path.join(paths["atlas_dir"], "nextbrain_segmentation", "assignments_nextbrain.csv")
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")
    #


    FF: As far as I understand, we dont use -scale_invnodevol nor -scale_invlenght,
    because f.e. -scale_invnodevol tells tck2connectome to scale each connection by
    the inverse of the parcel volume, but we don't need that because we have the
    external streamline weights from SIFT2 now
    Are these combinable?

    # It would be interesting to consider scaling by parcel volume or connection length with
    # -scale_invnodevol or -scale_invlength, respectively
    run_cmd([
        "tck2connectome",
        "-tck_weights_in", sift2_output,
        tckgen_output,
        parcels_coreg,
        connectome_csv,
        "-out_assignment", assignments_csv,
        "-symmetric", "-zero_diagonal",
        "-nthreads", str(102),
        "-force"
    ])

"""



if __name__ == '__main__':
    subject_dirs = get_subject_dirs("/media/nas/nikita/5p_connectome")
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        subject_id = os.path.basename(subj_dir)
        freesurfer_atlas_generation(paths, 60, subject_id)

