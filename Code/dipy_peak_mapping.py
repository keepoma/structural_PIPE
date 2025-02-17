import os
import nibabel as nib
import numpy as np
from nibabel.streamlines import load as load_streamlines
from dipy.tracking.streamline import values_from_volume
import matplotlib.pyplot as plt


def visualize_peak_components(peaks_3d, title="Peak Components", xlabel="Node Index", ylabel="Component Value"):
    """
    Visualizes the x, y, z components of the peak vector at each node.
    """

    n_nodes = peaks_3d.shape[0]
    x = np.arange(1, n_nodes + 1)

    plt.figure(figsize=(10, 6))
    plt.plot(x, peaks_3d[:, 0], label="x-component", color="red")
    plt.plot(x, peaks_3d[:, 1], label="y-component", color="green")
    plt.plot(x, peaks_3d[:, 2], label="z-component", color="blue")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()


def visualize_peak_magnitudes(peaks_3d, title="Peak Magnitudes", xlabel="Node Index", ylabel="Magnitude"):
    """
    Visualizes the magnitude (norm) of each peak vector along the streamline.
    """

    magnitudes = np.linalg.norm(peaks_3d, axis=1)
    x = np.arange(1, peaks_3d.shape[0] + 1)

    plt.figure(figsize=(10, 6))
    plt.plot(x, magnitudes, label="Peak Magnitude", color="purple")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()


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


# Visualize magnitude for the first streamline
first_streamline = np.array(sampled_peaks[0])
visualize_peak_components(first_streamline, title="Primary Peaks Along the First Streamline")
visualize_peak_magnitudes(first_streamline, title="Peak Magnitudes Along the First Streamline")

