# Generating peaks from previous data
sh2peaks wmfod_norm.mif fod_peaks.nii.gz -force

# Tract Segmentation, Bundle Start & End Regions, Tract Orientation Maps
export OMP_NUM_THREADS=1
TractSeg -i /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks.nii.gz \
         -o /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/tractseg_output \
         --nr_cpus 1
TractSeg -i /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks.nii.gz \
         -o /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/tractseg_output \
         --output_type endings_segmentation
TractSeg -i /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks.nii.gz \
         -o /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/tractseg_output \
         --output_type TOM
Tracking -i /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks.nii.gz \
         -o /Users/nikitakaruzin/MRI/projects/BATMAN/DWI/tractseg_output --nr_fibers 2000
