#!/bin/bash


# extract mean b0 from dMRI
for_each sample_598/* : dwiextract IN/5_dwi/dwi_den_unr_pre_unbia.mif -bzero IN/5_dwi/b0.nii.gz -fo;
for_each sample_598/* : mrmath IN/5_dwi/b0.nii.gz mean IN/2_nifti/mean_b0.nii.gz -axis 3 -fo;
for_each sample_598/* : rm IN/5_dwi/b0.nii.gz;
# mrinfo /media/nas/nikita/sample_598/raw/2_nifti/mean_b0.nii.gz

# 5tt from raw T1
for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/006_T1w_MPR IN/6_mif/t1_raw.mif -force;
for_each sample_598/* : 5ttgen fsl IN/6_mif/t1_raw.mif IN/6_mif/5tt_nocoreg -nthreads 50 -force;
for_each sample_598/* : mrconvert IN/6_mif/5tt_nocoreg.mif IN/2_nifti/5tt_nocoreg.nii.gz;
mrconvert T1_raw.mif T1_raw.nii.gz -nthreads 6 -force;

# Anatomically Constrained Tractography
for_each sample_598/* : flirt -in IN/2_nifti/mean_b0.nii.gz -ref IN/2_nifti/5tt_nocoreg.nii.gz -dof 6 -omat IN/7_mat/diff2struct_fsl.mat;
for_each sample_598/* : transformconvert IN/7_mat/diff2struct_fsl.mat IN/2_nifti/mean_b0.nii.gz IN/2_nifti/5tt_nocoreg.nii.gz flirt_import IN/7_mat/diff2struct_mrtrix.txt -nthreads 50 -force;
for_each sample_598/* : mrtransform IN/6_mif/5tt_nocoreg.mif -linear IN/7_mat/diff2struct_mrtrix.txt -inverse IN/6_mif/5tt_coreg.mif -nthreads 50 -force;
# check: mrview /media/nas/nikita/sample_598/raw/5_dwi/dwi_den_unr_pre_unbia.mif -overlay.load /media/nas/nikita/sample_598/raw/6_mif/t1_raw.mif - overlay.colormap 2 -overlay.load /media/nas/nikita/sample_598/raw/6_mif/t1_coreg.mif -overlay.colourmap 1

# Registering diffusion to structural space
flirt -in mean_b0_preprocessed.nii.gz -ref T1_raw.nii.gz -dof 6 -omat diff2struct_fsl.mat;
transformconvert diff2struct_fsl.mat mean_b0_preprocessed.nii.gz T1_raw.nii.gz flirt_import diff2struct_mrtrix.txt -nthreads 6 -force;
mrtransform T1_raw.mif -linear diff2struct_mrtrix.txt -inverse T1_coreg_diff2struct.mif -nthreads 6 -force;
mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg_diff2struct.mif -nthreads 6 -force;
# mrview dwi_den_unr_preproc_unbiased.mif -overlay.load T1_raw.mif -overlay.colourmap 2 -overlay.load T1_coreg_diff2struct.mif -overlay.colourmap 1

# Registering structural to diffusion space
mrconvert -strides 1,2,3,4 mean_b0_preprocessed.nii.gz mean_b0_preprocessed_newor.nii.gz -force;
mrconvert -strides 1,2,3,4 mean_b0_preprocessed.mif mean_b0_preprocessed_newor.mif -force;
mrconvert -strides 1,2,3 T1_raw.nii.gz T1_raw_newor.nii.gz -force;
mrconvert -strides 1,2,3 T1_raw.mif T1_raw_newor.mif -force;
mrconvert -strides 1,2,3,4 dwi_den_unr_preproc_unbiased.mif dwi_den_unr_preproc_unbiased_newor.mif -force;
flirt -in T1_raw_newor.nii.gz -ref mean_b0_preprocessed_newor.nii.gz -cost normmi -dof 6 -omat struct2diff_fsl.mat;
transformconvert struct2diff_fsl.mat T1_raw_newor.nii.gz mean_b0_preprocessed_newor.nii.gz flirt_import struct2diff_mrtrix.txt -nthreads 6 -force;
mrtransform T1_raw_newor.mif -linear struct2diff_mrtrix.txt T1_coreg_struct2diff.mif -nthreads 6 -force;
mrtransform 5tt_nocoreg.mif -linear struct2diff_mrtrix.txt 5tt_coreg_struct2diff.mif -nthreads 6 -force;
# mrview dwi_den_unr_preproc_unbiased_newor.mif -overlay.load T1_coreg_struct2diff.mif -overlay.colourmap 1 -overlay.opacity 0.5
# mrview dwi_den_unr_preproc_unbiased_newor.mif -overlay.load 5tt_coreg_struct2diff.mif -overlay.colourmap 1 -overlay.opacity 0.5


# Mask of streamline seeding
for_each sample_598/* : 5tt2gmwmi IN/6_mif/5tt_coreg.mif IN/6_mif/gmwmSeed_coreg.mif -fo;

# Tracks and SIFT
for_each sample_598/* : tckgen -act IN/6_mif/5tt_coreg.mif -backtrack -seed_gmwmi IN/6_mif/gmwmSeed_coreg.mif -select 10000000 IN/5_dwi/wm_norm.mif IN/8_tck/tracks_10mio.tck -nthreads 50 -fo;
for_each sample_598/* : tckedit IN/8_tck/tracks_10mio.tck -number 1000k IN/8_tck/smallerTracks_1000k.tck -nthreads 50 -fo;
# mrview '/media/nas/nikita/sample_598/raw/5_dwi/dwi_den_unr_pre_unbia.mif' -tractography.load '/media/nas/nikita/sample_598/raw/8_tck/smallerTracks_1000k.tck'
for_each sample_598/* : tcksift -act IN/6_mif/5tt_coreg.mif -term_number 1000000 IN/8_tck/tracks_10mio.tck IN/5_dwi/wm_norm.mif IN/8_tck/sift_1mio.tck -nthreads 50 -fo;
for_each sample_598/* : tckedit IN/8_tck/sift_1mio.tck -number 1000k IN/8_tck/smallerSIFT_1000k.tck -nthreads 50 -fo;
# mrview '/media/nas/nikita/sample_598/raw/5_dwi/dwi_den_unr_pre_unbia.mif' -tractography.load '/media/nas/nikita/sample_598/raw/8_tck/smallerSIFT_1000k.tck'

# Region of interest localization
for_each sample_598/* : tckedit -include -26.5,33.95,27.41,3 IN/8_tck/sift_1mio.tck IN/8_tck/FA.tck -nthreads 50 -fo;
for_each sample_598/* : mrtransform IN/6_mif/T1_raw.mif -linear IN/7_mat/diff2struct_mrtrix.txt -inverse IN/6_mif/t1_coreg.mif -nthreads 50 -force;
# mrview '/media/nas/nikita/sample_598/raw/6_mif/t1_coreg.mif' -tractography.load '/media/nas/nikita/sample_598/raw/8_tck/FA.tck'

# Freesurfer/HCP-based atlas generation
recon-all -s subject_598 -i /media/nas/nikita/sample_598/raw/2_nifti/T1.nii.gz -all -threads $(nproc --ignore=10)
# at group level consider: recon-all -s <subjectName> -qcache
# Mapping of HCP MMP 1.0 atlas from fsaverage onto subject
mri_surf2surf --srcsubject fsaverage --trgsubject subject_598 --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.hcpmmp1.annot --tval $SUBJECTS_DIR/subject_598/label/lh.hcpmmp1.annot;
mri_surf2surf --srcsubject fsaverage --trgsubject subject_598 --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.hcpmmp1.annot --tval $SUBJECTS_DIR/subject_598/label/rh.hcpmmp1.annot;
# Mapping HCP annotations onto volumetric image, subcortical segmentation.
mri_aparc2aseg --old-ribbon --s subject_598 --annot hcpmmp1 --o /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.mgz;
mrconvert -datatype uint32 /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.mgz /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.mif;
# Integer ordering
labelconvert /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.mif /home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_original.txt /home/nikita/anaconda3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1_parcels_nocoreg.mif;
# Registering atlas-based volumetric parcellation to dMRI
mrtransform /home/nikita/Structural_PIPE/hcpmmp1_parcels_nocoreg.mif -linear /media/nas/nikita/sample_598/raw/7_mat/diff2struct_mrtrix.txt -inverse -datatype uint32 /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1_parcels_coreg.mif
#mrview '/media/nas/nikita/sample_598/raw/6_mif/t1_coreg.mif' -overlay.load '/media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1_parcels_coreg.mif' -overlay.colourmap random -overlay.opacity 0.4

# Matrix Generation
tck2connectome -symmetric -zero_diagonal -scale_invnodevol /media/nas/nikita/sample_598/raw/8_tck/sift_1mio.tck /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1_parcels_coreg.mif /media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.csv -out_assignment /media/nas/nikita/sample_598/raw/9_atlas/assignments_hcpmmp1.csv

# Extrating streamlines between 15 atlas regions
# I localized 15 streamlines of interest for the speech regions of the left hemisphere and generates
# 100+ streamlines to analyze their connections. Script in "custom_speechR_from_global.py"



