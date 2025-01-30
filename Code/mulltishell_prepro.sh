#!/bin/bash

for_each * : mrconvert -strides 1,2,3 IN/1_raw/006_T1w_MPR IN/2_nifti/t1.nii.gz -force;

for_each * : mrconvert -strides 1,2,3 IN/1_raw/008_T2w_SPC IN/2_nifti/t2.nii.gz -force;

for_each * : mrconvert -strides 1,2,3 IN/1_raw/009_t2_space_dark-fluid_sag_p2_iso_0_8 IN/2_nifti/t2_df.nii.gz -force;

for_each * : mrconvert -strides 1,2,3,4 IN/1_raw/016_dMRI_dir98_AP IN/5_dwi/dwi_ap.mif -force;

for_each * : mrconvert -strides 1,2,3,4 IN/1_raw/019_dMRI_dir98_PA IN/5_dwi/dwi_pa.mif -force;


#preprocess
#multi-shell with ap/pa:

for_each * : mrcat IN/5_dwi/dwi_ap.mif IN/5_dwi/dwi_pa.mif IN/5_dwi/dwi_all.mif -axis 3;
for_each * : dwidenoise IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -noise IN/5_dwi/noise.mif;
for_each * : mrcalc IN/5_dwi/dwi_all.mif IN/5_dwi/dwi_den.mif -subtract IN/5_dwi/residual.mif;

#Unringing
#The “axes” option must be adjusted to your dataset: With this option, you inform the algorithm of the plane in which you acquired your data: –axes 0,1 means you acquired axial slices; -axes 0,2 refers to coronal slices and –axes 1,2 to sagittal slices!
for_each * : mrdegibbs IN/5_dwi/dwi_den.mif IN/5_dwi/dwi_den_unr.mif -axes 0,1;

#motion and distortion correction
for_each * : dwifslpreproc -rpe_all -pe_dir AP IN/5_dwi/dwi_den_unr.mif IN/5_dwi/dwi_den_unr_pre.mif;

#Bias field correction
for_each * : dwibiascorrect ants IN/5_dwi/dwi_den_unr_pre.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif -bias IN/5_dwi/bias.mif;

#generate mask
for_each * : dwi2mask IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif;

#skull stripping dwi
for_each * : mrcalc IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/mask.mif -mult IN/5_dwi/dwi_den_unr_pre_unbia_skull.mif

#Fiber orientation distribution
for_each * : dwi2response dhollander IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/gm.txt IN/5_dwi/csf.txt;
for_each * : dwi2fod msmt_csd -mask IN/5_dwi/mask.mif IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/wm.txt IN/5_dwi/wm.mif IN/5_dwi/gm.txt IN/5_dwi/gm.mif IN/5_dwi/csf.txt IN/5_dwi/csf.mif;

# extract mean b0 from dMRI
for_each * : dwiextract IN/5_dwi/dwi_den_unr_pre_unbia.mif -bzero IN/5_dwi/b0.nii.gz -fo;
for_each * : mrmath IN/5_dwi/b0.nii.gz mean IN/2_nifti/b0.nii.gz -axis 3 -fo;
for_each * : rm IN/5_dwi/b0.nii.gz;
