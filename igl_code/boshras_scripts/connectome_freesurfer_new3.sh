#!/bin/bash

#construct the connectome (currently 40M, might change to 20M - 60M streamlines)
#for_each subjects/1_patients_3/* : tckgen IN/5_DWI/wm.mif IN/9_connectome_fs/40M.tck -act IN/8_5tt_and_lesion/5tt_hsvs_vbg_lesion.mif -backtrack -crop_at_gmwmi -seed_dynamic IN/5_DWI/wm.mif -maxlength 250 -select 40M -power 0.33

#sift the tractogram
#Use SIFT2 to determine streamline weights
for_each subjects/1_patients_recheck/* : tcksift2 IN/9_connectome_fs/40M_lesioned.tck IN/5_DWI/wm.mif IN/9_connectome_fs/weights_lesioned.csv -act IN/8_5tt_and_lesion/5tt_hsvs_vbg_lesion.mif -out_mu IN/9_connectome_fs/mu_lesioned.txt -out_coeffs IN/9_connectome_fs/tck_coeffs_lesioned.txt

# Generate TDIs
for_each subjects/1_patients_recheck/* : tckmap IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/tdi_hires_lesioned.mif -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv -vox 0.25 -datatype uint16

for_each subjects/1_patients_recheck/* : tckmap IN/9_connectome_fs/40M_lesioned.tck - -template IN/9_connectome_fs/brain.mgz -precise | mrcalc - $(cat IN/9_connectome_fs/mu_lesioned.txt) -mult IN/9_connectome_fs/tdi_T1_lesioned.mif

#convert labels for parcellation
#for_each subjects/1_patients_recheck/* : labelconvert IN/9_connectome_fs/aparc+aseg.mgz /usr/local/freesurfer/FreeSurferColorLUT.txt /home/adt/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt IN/9_connectome_fs/parc_fs.nii.gz

#for_each subjects/1_patients_recheck/* : label2mesh IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/nodes.obj 
#for_each subjects/1_patients_recheck/* : meshfilter IN/9_connectome_fs/nodes.obj smooth IN/9_connectome_fs/nodes_smooth.obj

#construct connectome
for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/connectome_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv -out_assignments IN/9_connectome_fs/assignments_lesioned.txt -symmetric -zero_diagonal


#generate data for visualisation in mrtrix

for_each subjects/1_patients_recheck/* : connectome2tck IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/assignments_lesioned.txt IN/9_connectome_fs/exemplars_lesioned.tck -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv -exemplars IN/9_connectome_fs/parc_fs.nii.gz -files single 

##calculate tensors
#for_each subjects/1_patients_recheck/* : dwi2tensor IN/5_DWI/dwi_den_unr_pre_unbia_up.mif -mask IN/5_DWI/mask_up.mif IN/5_DWI/tensors.mif

## calculate dMRI metrics
#for_each subjects/1_patients_recheck/* : tensor2metric IN/5_DWI/tensors.mif -mask IN/5_DWI/mask_up.mif -fa IN/5_DWI/fa.mif
#for_each subjects/1_patients_recheck/* : tensor2metric IN/5_DWI/tensors.mif -mask IN/5_DWI/mask_up.mif -adc IN/5_DWI/adc.mif
#for_each subjects/1_patients_recheck/* : tensor2metric IN/5_DWI/tensors.mif -mask IN/5_DWI/mask_up.mif -ad IN/5_DWI/ad.mif
#for_each subjects/1_patients_recheck/* : tensor2metric IN/5_DWI/tensors.mif -mask IN/5_DWI/mask_up.mif -rd IN/5_DWI/rd.mif

## connectome edge: mean dMRI metrics along tracks 
for_each subjects/1_patients_recheck/* : tcksample IN/9_connectome_fs/40M_lesioned.tck IN/5_DWI/fa.mif IN/9_connectome_fs/tck_meanFA_lesioned.csv -stat_tck mean 
for_each subjects/1_patients_recheck/* : tcksample IN/9_connectome_fs/40M_lesioned.tck IN/5_DWI/adc.mif IN/9_connectome_fs/tck_meanADC_lesioned.csv -stat_tck mean 
for_each subjects/1_patients_recheck/* : tcksample IN/9_connectome_fs/40M_lesioned.tck IN/5_DWI/ad.mif IN/9_connectome_fs/tck_meanAD_lesioned.csv -stat_tck mean
for_each subjects/1_patients_recheck/* : tcksample IN/9_connectome_fs/40M_lesioned.tck IN/5_DWI/rd.mif IN/9_connectome_fs/tck_meanRD_lesioned.csv -stat_tck mean 

# dMRI metrics connectome
# Here, a connectome matrix that is "weighted by FA" is generated in multiple steps: firstly, for each streamline, the value of the underlying FA image is sampled at each vertex, and the mean of these values is calculated to produce a single scalar value of "mean FA" per streamline; then, as each streamline is assigned to nodes within the connectome, the magnitude of the contribution of that streamline to the matrix is multiplied by the mean FA value calculated prior for that streamline; finally, for each connectome edge, across the values of "mean FA" that were contributed by all of the streamlines assigned to that particular edge, the mean value is calculated

for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/connectome_fa_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv  -scale_file IN/9_connectome_fs/tck_meanFA_lesioned.csv -stat_edge mean -symmetric -zero_diagonal 

for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/connectome_adc_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv  -scale_file IN/9_connectome_fs/tck_meanADC_lesioned.csv -stat_edge mean -symmetric -zero_diagonal 

for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/connectome_ad_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv  -scale_file IN/9_connectome_fs/tck_meanAD_lesioned.csv -stat_edge mean -symmetric -zero_diagonal 

for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/connectome_rd_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv  -scale_file IN/9_connectome_fs/tck_meanRD_lesioned.csv -stat_edge mean -symmetric -zero_diagonal 

#mean edge length (Generate a matrix consisting of the mean streamline length between each node pair)
for_each subjects/1_patients_recheck/* : tck2connectome IN/9_connectome_fs/40M_lesioned.tck IN/9_connectome_fs/parc_fs.nii.gz IN/9_connectome_fs/meanlength_lesioned.csv -tck_weights_in IN/9_connectome_fs/weights_lesioned.csv -scale_length -stat_edge mean -symmetric -zero_diagonal



#single
#tckgen 5_DWI/wm.mif 9_connectome_fs/40M.tck -act 8_5tt_and_lesion/5tt_hsvs_vbg_lesion.mif -backtrack -crop_at_gmwmi -seed_dynamic 5_DWI/wm.mif -maxlength 250 -select 40M -power 0.33

#tcksift2 9_connectome_fs/40M.tck 5_DWI/wm.mif 9_connectome_fs/weights.csv -act 8_5tt_and_lesion/5tt_hsvs_vbg_lesion.mif -out_mu 9_connectome_fs/mu.txt -out_coeffs 9_connectome_fs/tck_coeffs.txt

#tckmap 9_connectome_fs/40M.tck 9_connectome_fs/tdi_hires.mif -tck_weights_in 9_connectome_fs/weights.csv -vox 0.25 -datatype uint16

#tckmap 9_connectome_fs/40M.tck - -template 9_connectome_fs/brain.mgz -precise | mrcalc - $(cat 9_connectome_fs/mu.txt) -mult 9_connectome_fs/tdi_T1.mif

#labelconvert 9_connectome_fs/aparc+aseg.mgz /usr/local/freesurfer/FreeSurferColorLUT.txt  /home/adt/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt 9_connectome_fs/parc_fs.nii.gz

#label2mesh 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/nodes.obj 

#meshfilter 9_connectome_fs/nodes.obj smooth 9_connectome_fs/nodes_smooth.obj

#tck2connectome 9_connectome_fs/40M.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/connectome.csv -tck_weights_in 9_connectome_fs/weights.csv -out_assignments 9_connectome_fs/assignments.txt -symmetric -zero_diagonal

#connectome2tck 9_connectome_fs/40M.tck 9_connectome_fs/assignments.txt 9_connectome_fs/exemplars.tck -tck_weights_in 9_connectome_fs/weights.csv -exemplars 9_connectome_fs/parc_fs.nii.gz -files single

#tcksample 9_connectome_fs/40M.tck 5_DWI/fa.mif 9_connectome_fs/tck_meanFA.csv -stat_tck mean 
#tcksample 9_connectome_fs/40M.tck 5_DWI/adc.mif 9_connectome_fs/tck_meanADC.csv -stat_tck mean 
#tcksample 9_connectome_fs/40M.tck 5_DWI/ad.mif 9_connectome_fs/tck_meanAD.csv -stat_tck mean
#tcksample 9_connectome_fs/40M.tck 5_DWI/rd.mif 9_connectome_fs/tck_meanRD.csv -stat_tck mean 



#tck2connectome 9_connectome_fs/40M.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/connectome_fa.csv -tck_weights_in 9_connectome_fs/weights.csv  -scale_file 9_connectome_fs/tck_meanFA.csv -stat_edge mean -symmetric -zero_diagonal -fo
#tck2connectome 9_connectome_fs/40M.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/connectome_adc.csv -tck_weights_in 9_connectome_fs/weights.csv  -scale_file 9_connectome_fs/tck_meanADC.csv -stat_edge mean -symmetric -zero_diagonal -fo
#tck2connectome 9_connectome_fs/40M.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/connectome_ad.csv -tck_weights_in 9_connectome_fs/weights.csv  -scale_file 9_connectome_fs/tck_meanAD.csv -stat_edge mean -symmetric -zero_diagonal -fo
#tck2connectome 9_connectome_fs/40M.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/connectome_rd.csv -tck_weights_in 9_connectome_fs/weights.csv  -scale_file 9_connectome_fs/tck_meanRD.csv -stat_edge mean -symmetric -zero_diagonal -fo

#mean edge length (Generate a matrix consisting of the mean streamline length between each node pair)
#tck2connectome 9_connectome_fs/40M_lesioned.tck 9_connectome_fs/parc_fs.nii.gz 9_connectome_fs/meanlength.csv -tck_weights_in 9_connectome_fs/weights.csv -scale_length -stat_edge mean -symmetric -zero_diagonal -fo






