#!/bin/bash
#FOD based metrics

#Perform segmentation of continuous Fibre Orientation Distributions (FODs) to produce discrete fixels with use of option flag '-afd image'

for_each subjects/* : fod2fixel -afd fd.mif IN/5_dwi/wm.mif IN/5_dwi/fd_output -fo


#Compute fixel map from a tractogram,left hemisphere

for_each subjects/* : tck2fixel IN/9_connectome_fs/40M.tck IN/5_dwi/fd_output IN/5_dwi/fd_output tck2fixel_input_mean_fd.mif


#Generate a fixel mask, remove all fixels where tck2fixel.mif==0

for_each subjects/* : fixelcrop IN/5_dwi/fd_output IN/5_dwi/fd_output/tck2fixel_input_mean_fd.mif IN/5_dwi/mean_fd


#Convert a fixel-based sparse-data image into a scalar image

for_each subjects/* : fixel2voxel IN/5_dwi/mean_fd/fd.mif mean IN/5_dwi/mean_fd/fixel2voxel_input_fd_weighted.mif -weighted IN/5_dwi/mean_fd/fd.mif
