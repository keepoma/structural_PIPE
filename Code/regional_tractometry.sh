# Generating peaks from previous data
sh2peaks wmfod_norm.mif fod_peaks.nii.gz -force

# Tract Segmentation, Bundle Start & End Regions, Tract Orientation Maps
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type tract_segmentation
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type endings_segmentation
TractSeg -i '/media/nas/nikita/sample_598/raw/2_nifti/fod_peaks.nii.gz' --output_type TOM

