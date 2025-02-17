import os
import nibabel as nib
import numpy as np
from nibabel.streamlines import load as load_streamlines
from dipy.tracking.streamline import values_from_volume


# Paths
peaks_path_group_RF = "/media/nas/nikita/test_study2_1sub/test_302/raw/2_nifti/fod_peaks_group_RF.nii.gz"
AF_left_resampled = "/media/nas/nikita/test_study2_1sub/test_302/tractseg_output/FOD_iFOD2_trackings/AF_left_N100.tck"

# Loading peaks image
peaks_img = nib.load(peaks_path_group_RF)
peaks_data = peaks_img.get_fdata()
affine = peaks_img.affine
print("Peaks image shape:", peaks_data.shape)

# Extract the first three components.
primary_peaks = peaks_data[..., :3]

# Load the streamlines from the TCK file
sft = load_streamlines(AF_left_resampled)
streamlines = list(sft.streamlines)

# Map the primary peaks along each resampled streamline
sampled_peaks = values_from_volume(primary_peaks, streamlines, affine)

# Compute the average primary peak along the first streamline:
avg_peak_first = np.mean(sampled_peaks[0], axis=0)
print("Average primary peak for first streamline:", avg_peak_first)
