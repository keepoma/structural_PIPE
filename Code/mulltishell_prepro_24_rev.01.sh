#!/bin/bash

#with "for_each" "-nthreads #" will have different results if used before or after the colon
#for_each study/* -nthreads 2 : dwidostuff IN/dwi.mif IN/dwi_stuffdone.mif -nthreads 4, This would instruct for_each to always have two jobs running in parallel, each of which will be explicitly instructed to use four threads.

for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/006_T1w_MPR IN/2_nifti/t1.nii.gz -force;

for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/008_T2w_SPC IN/2_nifti/t2.nii.gz -force;

for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/009_t2_space_dark-fluid_sag_p2_iso_0_8 IN/2_nifti/t2_df.nii.gz -force;

for_each sample_598/* : mrconvert -strides 1,2,3,4 IN/1_raw/016_dMRI_dir98_AP IN/5_dwi/dwi_ap.mif -force;

for_each sample_598/* : mrconvert -strides 1,2,3,4 IN/1_raw/019_dMRI_dir98_PA IN/5_dwi/dwi_pa.mif -force;


#preprocess
#multi-shell with ap/pa:

for_each sample_598/* : mrcat IN/5_dwi/dwi_ap.mif IN/5_dwi/dwi_pa.mif IN/5_dwi/dwi_all.mif -axis 3;
for_each sample_598/* : dwidenoise IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -noise IN/5_dwi/noise.mif;
for_each sample_598/* : mrcalc IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -subtract IN/5_dwi/residual.mif;

#Unringing
#The “axes” option must be adjusted to your dataset: With this option, you inform the algorithm of the plane in which you acquired your data: –axes 0,1 means you acquired axial slices; -axes 0,2 refers to coronal slices and –axes 1,2 to sagittal slices!
for_each sample_598/* : mrdegibbs IN/5_dwi/dwi_den.mif IN/5_dwi/dwi_den_unr.mif -axes 0,1;

#motion and distortion correction
for_each sample_598/* : dwifslpreproc -rpe_all -pe_dir AP IN/5_dwi/dwi_den_unr.mif IN/5_dwi/dwi_den_unr_pre.mif;

#Bias field correction
for_each sample_598/* : dwibiascorrect ants IN/5_dwi/dwi_den_unr_pre.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif -bias IN/5_dwi/bias.mif;

#generate mask
for_each sample_598/* : dwi2mask IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif;

#skull stripping dwi
for_each sample_598/* : mrcalc IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif -mult IN/5_dwi/dwi_den_unr_pre_unbia_skull.mif;

#Fiber orientation distribution
for_each sample_598/* : dwi2response dhollander IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/gm.txt IN/5_dwi/csf.txt;
for_each sample_598/* : dwi2fod msmt_csd -mask IN/5_dwi/mask.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/wm.mif IN/5_dwi/gm.txt IN/5_dwi/gm.mif IN/5_dwi/csf.txt IN/5_dwi/csf.mif;
#check: using wm FOD
# mrconvert -coord 3 0 wm.mif - | mrcat csf.mif gm.mif - vf.mif mrview vf.mif -odf.load_sh wm.mif

#Intensity Normalization
for_each sample_598/* : mtnormalise IN/5_dwi/wm.mif IN/5_dwi/wm_norm.mif IN/5_dwi/gm.mif IN/5_dwi/gm_norm.mif IN/5_dwi/csf.mif IN/5_dwi/csf_norm.mif -mask IN/5_dwi/mask.mif;

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
# https://www.youtube.com/watch?v=A2ZyGE5BcfE&list=PLIQIswOrUH68Zi9SVDAdcUExpq2i6A2eD&index=7 - revisar
# THESE 2 WORKED fslreorient2std t1.nii.gz t1_std.nii.gz
# fslreorient2std mean_b0.nii.gz mean_b0_std.nii.gz
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
recon-all -s subject_598 -i /media/nas/nikita/sample_598/raw/2_nifti/T1.nii.gz -all

