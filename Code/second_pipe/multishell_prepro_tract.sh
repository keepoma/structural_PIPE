#!/bin/bash


#with "for_each" "-nthreads #" will have different results if used before or after the colon
#for_each study/* -nthreads 2 : dwidostuff IN/dwi.mif IN/dwi_stuffdone.mif -nthreads 4, This would instruct for_each to always have two jobs running in parallel, each of which will be explicitly instructed to use four threads.

for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/006_T1w_MPR IN/2_nifti/t1.nii.gz -force;
for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/008_T2w_SPC IN/2_nifti/t2.nii.gz -force;
for_each sample_598/* : mrconvert -strides 1,2,3 IN/1_raw/009_t2_space_dark-fluid_sag_p2_iso_0_8 IN/2_nifti/t2_df.nii.gz -force;
for_each sample_598/* : mrconvert -strides 1,2,3,4 IN/1_raw/016_dMRI_dir98_AP IN/5_dwi/dwi_ap.mif -force;
for_each sample_598/* : mrconvert -strides 1,2,3,4 IN/1_raw/019_dMRI_dir98_PA IN/5_dwi/dwi_pa.mif -force;

# Preprocess
# Multi-shell with ap/pa:
for_each sample_598/* : mrcat IN/5_dwi/dwi_ap.mif IN/5_dwi/dwi_pa.mif IN/5_dwi/dwi_all.mif -axis 3;
for_each sample_598/* : dwidenoise IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -noise IN/5_dwi/noise.mif;
for_each sample_598/* : mrcalc IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -subtract IN/5_dwi/residual.mif;

# Unringing
# The “axes” option must be adjusted to your dataset: With this option, you inform the algorithm of the plane in which you acquired your data: –axes 0,1 means you acquired axial slices; -axes 0,2 refers to coronal slices and –axes 1,2 to sagittal slices!
for_each sample_598/* : mrdegibbs IN/5_dwi/dwi_den.mif IN/5_dwi/dwi_den_unr.mif -axes 0,1;

# Motion and distortion correction
for_each sample_598/* : dwifslpreproc -rpe_all -pe_dir AP IN/5_dwi/dwi_den_unr.mif IN/5_dwi/dwi_den_unr_pre.mif;

# Bias field correction
for_each sample_598/* : dwibiascorrect ants IN/5_dwi/dwi_den_unr_pre.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif -bias IN/5_dwi/bias.mif;

# Generate mask
for_each sample_598/* : dwi2mask IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif;

# Skull stripping dwi
for_each sample_598/* : mrcalc IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif -mult IN/5_dwi/dwi_den_unr_pre_unbia_skull.mif;

# Response Function estimation and fiber orientation distribution
for_each sample_598/* : dwi2response dhollander IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/gm.txt IN/5_dwi/csf.txt;
for_each sample_598/* : dwi2fod msmt_csd -mask IN/5_dwi/mask.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/wm.mif IN/5_dwi/gm.txt IN/5_dwi/gm.mif IN/5_dwi/csf.txt IN/5_dwi/csf.mif;
#check: using wm FOD
# mrconvert -coord 3 0 wm.mif - | mrcat csf.mif gm.mif - vf.mif mrview vf.mif -odf.load_sh wm.mif

# Intensity Normalization
for_each sample_598/* : mtnormalise IN/5_dwi/wm.mif IN/5_dwi/wm_norm.mif IN/5_dwi/gm.mif IN/5_dwi/gm_norm.mif IN/5_dwi/csf.mif IN/5_dwi/csf_norm.mif -mask IN/5_dwi/mask.mif;








