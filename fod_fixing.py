import nibabel as nib
import numpy as np

# Load file
peaks_nii = nib.load("/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks.nii.gz")
peaks_data = peaks_nii.get_fdata()

# Extract first 3 components (X, Y, Z)
peaks_data_fixed = peaks_data[:, :, :, :3]

# Convert to float32
peaks_data_fixed = peaks_data_fixed.astype(np.float32)

# Save as new NIfTI file
fixed_nii = nib.Nifti1Image(peaks_data_fixed, affine=peaks_nii.affine, header=peaks_nii.header)
nib.save(fixed_nii, "/Users/nikitakaruzin/MRI/projects/BATMAN/DWI/fod_peaks_fixed.nii.gz")