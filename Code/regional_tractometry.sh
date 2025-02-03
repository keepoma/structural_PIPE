# DTI
for_each sample_598/* : dwi2tensor IN/5_dwi/dwi_den_unr_pre_unbia.mif IN/5_dwi/tensor.mif;
for_each sample_598/* : tensor2metric IN/5_dwi/tensor.mif -adc IN/5_dwi/adc.mif;
for_each sample_598/* : tensor2metric IN/5_dwi/tensor.mif -fa IN/5_dwi/fa.mif;


# Generating peaks from previous data
sh2peaks wmfod_norm.mif fod_peaks.nii.gz -force

# Tract Segmentation, Bundle Start & End Regions, Tract Orientation Maps
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type tract_segmentation
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type endings_segmentation
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type TOM

