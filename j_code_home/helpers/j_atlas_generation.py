import os
import re
from helpers.j_helpers import run_cmd, get_subject_paths, get_subject_dirs


"""
Functions for atlas generation.
Currently based on:
    -NextBrain
    -Freesurfer
"""


def freesurfer_atlas_generation(paths, nthreads, subject_id_ses):
    """
    Generate the Freesurfer/HCP-based atlas and perform label conversion.
    """
    atlas_dir = paths["atlas_dir"]
    parcels_coreg = os.path.join(atlas_dir, "hcpmmp1_parcels_coreg.mif")
    subjects_dir = "/home/nikita/freesurfer/subjects"
    subject_dir = os.path.join(subjects_dir, subject_id_ses)
    if not os.path.isdir(subject_dir):
        os.makedirs(paths["atlas_dir"], exist_ok=True)
        # Run Freesurfer's recon-all using the subject ID
        # use command rm -rf $SUBJECTS_DIR/subject_id if rerunning
        for filename in os.listdir(paths["anat_dir"]):
            if re.search(r'T1w?\.nii\.gz', filename):
                t1_nii = os.path.join(paths["anat_dir"], filename)
        run_cmd([
            "recon-all",
            "-s", subject_id_ses,
            "-i", t1_nii,
            "-all",
            "-threads", str(nthreads)
        ])
    else:
        print("Found existing directory in freesurfer, delete file if necessary, skipping step")

    # Map HCP MMP 1.0 atlas from fsaverage onto the subject
    subjects_dir = os.environ.get("SUBJECTS_DIR", "")
    lh_src = os.path.join(subjects_dir, "fsaverage", "label", "lh.hcpmmp1.annot")
    lh_trg = os.path.join(subjects_dir, subject_id_ses, "label", "lh.hcpmmp1.annot")
    rh_src = os.path.join(subjects_dir, "fsaverage", "label", "rh.hcpmmp1.annot")
    rh_trg = os.path.join(subjects_dir, subject_id_ses, "label", "rh.hcpmmp1.annot")

    if not os.path.exists(parcels_coreg):
        run_cmd([
            "mri_surf2surf",
            "--srcsubject", "fsaverage",
            "--trgsubject", subject_id_ses,
            "--hemi", "lh",
            "--sval-annot", lh_src,
            "--tval", lh_trg
        ])

        run_cmd([
            "mri_surf2surf",
            "--srcsubject", "fsaverage",
            "--trgsubject", subject_id_ses,
            "--hemi", "rh",
            "--sval-annot", rh_src,
            "--tval", rh_trg
        ])


        # Create atlas outputs in the designated atlas folder
        atlas_mgz = os.path.join(atlas_dir, "hcpmmp1.mgz")
        atlas_mif = os.path.join(atlas_dir, "hcpmmp1.mif")

        run_cmd([
            "mri_aparc2aseg",
            "--old-ribbon",
            "--s", subject_id_ses,
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
    else:
        print("Parcels coreg exists, skipping atlas")


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

