#!/bin/bash


# extract mean b0 from dMRI
for_each sample_598/* : dwiextract IN/5_dwi/dwi_den_unr_pre_unbia.mif -bzero IN/5_dwi/b0.nii.gz -fo;
for_each sample_598/* : mrmath IN/5_dwi/b0.nii.gz mean IN/2_nifti/mean_b0.nii.gz -axis 3 -fo;
for_each sample_598/* : rm IN/5_dwi/b0.nii.gz;
# mrinfo /media/nas/nikita/sample_598/raw/2_nifti/mean_b0.nii.gz

# 5tt T1
for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/006_T1w_MPR IN/6_mif/t1_raw.mif -force;
for_each sample_598/* : 5ttgen fsl IN/6_mif/t1_raw.mif IN/6_mif/5tt_nocoreg -nthreads 50 -force;
for_each sample_598/* : mrconvert IN/6_mif/5tt_nocoreg.mif IN/2_nifti/5tt_nocoreg.nii.gz;

# Anatomically Constrained Tractography
for_each sample_598/* : flirt -in IN/2_nifti/mean_b0.nii.gz -ref IN/2_nifti/5tt_nocoreg.nii.gz -dof 6 -omat IN/7_mat/diff2struct_fsl.mat;
for_each sample_598/* : transformconvert IN/7_mat/diff2struct_fsl.mat IN/2_nifti/mean_b0.nii.gz IN/2_nifti/5tt_nocoreg.nii.gz flirt_import IN/7_mat/diff2struct_mrtrix.txt -nthreads 50 -force;
for_each sample_598/* : mrtransform IN/6_mif/5tt_nocoreg.mif -linear IN/7_mat/diff2struct_mrtrix.txt -inverse IN/6_mif/5tt_coreg.mif -nthreads 50 -force;
# check: mrview /media/nas/nikita/sample_598/raw/5_dwi/dwi_den_unr_pre_unbia.mif -overlay.load /media/nas/nikita/sample_598/raw/6_mif/t1_raw.mif - overlay.colormap 2 -overlay.load /media/nas/nikita/sample_598/raw/6_mif/t1_coreg.mif -overlay.colourmap 1

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

#Freesurfer
recon-all -s subject_598 -i /media/nas/nikita/sample_598/raw/2_nifti/T1.nii.gz -all -threads $(nproc --ignore=10)