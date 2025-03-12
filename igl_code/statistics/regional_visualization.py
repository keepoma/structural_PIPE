import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def visualize_peak_length(peaks_txt):
    """
    Generates mean peak amplitude from 2000 streamlines across 100 points
    """

    # Load data as 1D array
    data_1d = np.loadtxt(peaks_txt)
    n_points = 100
    n_streamlines = 2000

    # Reshape data to (n_streamlines, n_points)
    data_2d = data_1d.reshape(n_streamlines, n_points)


    # Mean and SD plot
    mean_values = data_2d.mean(axis=0)
    std_values = data_2d.std(axis=0)

    plt.figure(figsize=(8, 5))
    plt.plot(mean_values, label="Mean")
    plt.fill_between(
        np.arange(n_points),
        mean_values - std_values,
        mean_values + std_values,
        alpha=0.2, label="±1 SD"
    )

    plt.title("Mean ± SD of Peak Values")
    plt.xlabel("Node index")
    plt.ylabel("Peak amplitude")
    plt.legend()
    plt.show()


def visualize_multi_peak_length_side_by_side(peaks_txt_files):
    """
    Visualizes multiple peak amplitude profiles in one figure with
    two subplots:
      1) Raw peak amplitude
      2) Absolute-value peak amplitude
    """

    plt.style.use('fast')
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

    colors = ["blue", "red", "green", "orange"]
    n_points = 100
    n_streamlines = 2000
    x = np.arange(n_points)

    # -- Subplot 1: Raw Peak Amplitudes --
    for idx, txt_file in enumerate(peaks_txt_files):
        data_1d = np.loadtxt(txt_file)
        data_2d = data_1d.reshape(n_streamlines, n_points)

        mean_values = data_2d.mean(axis=0)
        std_values  = data_2d.std(axis=0)

        color     = colors[idx % len(colors)]
        label_txt = os.path.basename(txt_file)

        axes[0].plot(x, mean_values,
                     color=color, linewidth=2,
                     label=f"{label_txt} Mean")
        axes[0].fill_between(x,
                             mean_values - std_values,
                             mean_values + std_values,
                             color=color, alpha=0.2,
                             label=f"{label_txt} ±1 SD")

    axes[0].set_title("Mean ± SD of Peak Values (Raw)")
    axes[0].set_xlabel("Node Index")
    axes[0].set_ylabel("Peak Amplitude")
    axes[0].legend(loc='lower left')
    axes[0].grid(True)

    # -- Subplot 2: Absolute Peak Amplitudes --
    for idx, txt_file in enumerate(peaks_txt_files):
        data_1d = np.loadtxt(txt_file)
        data_2d = data_1d.reshape(n_streamlines, n_points)

        # Take absolute value
        data_2d_abs = np.abs(data_2d)

        mean_values_abs = data_2d_abs.mean(axis=0)
        std_values_abs  = data_2d_abs.std(axis=0)

        color     = colors[idx % len(colors)]
        label_txt = os.path.basename(txt_file)

        axes[1].plot(x, mean_values_abs,
                     color=color, linewidth=2,
                     label=f"{label_txt} Mean (Abs)")
        axes[1].fill_between(x,
                             mean_values_abs - std_values_abs,
                             mean_values_abs + std_values_abs,
                             color=color, alpha=0.2,
                             label=f"{label_txt} ±1 SD (Abs)")

    axes[1].set_title("Mean ± SD of Peak Values (Absolute)")
    axes[1].set_xlabel("Node Index")
    axes[1].set_ylabel("Peak Amplitude")
    axes[1].legend(loc='lower left')
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()




def visualize_fa(csv_file):
    """
    Visualizes the FA profile along a tract from a CSV file.
    """

    # Read the FA CSV file into a df
    data = pd.read_csv(csv_file, skiprows=1, header=None, sep='\s+')

    plt.figure(figsize=(10, 5))

    mean_values = data.mean(axis=0)
    std_values = data.std(axis=0)
    x = range(1, len(mean_values) + 1)
    plt.plot(x, mean_values, label="Mean FA", color="blue", linewidth=2)
    plt.fill_between(x, mean_values - std_values, mean_values + std_values,
                     color="blue", alpha=0.3, label="Std. Deviation")

    plt.xlabel("Along-Tract Position (resampled)")
    plt.ylabel("FA Value")
    plt.title(os.path.basename(csv_file))
    plt.legend()
    plt.grid(True)
    plt.show()


def visualize_multi_fa(dirs, csv_files):
    """
    Visualizes multiple FA profiles in one plot
    """
    plt.style.use('fast')
    plt.figure(figsize=(12, 6))

    # Define colors for up to four files
    colors = ["blue", "red", "green", "orange"]

    for idx, csv_file in enumerate(csv_files):
        data = pd.read_csv(csv_file, skiprows=1, header=None, sep='\s+')
        mean_values = data.mean(axis=0)
        std_values = data.std(axis=0)
        x = range(1, len(mean_values) + 1)

        color = colors[idx % len(colors)]
        # Plot mean FA
        plt.plot(x, mean_values, label=f"{os.path.basename(csv_file)} Mean FA", color=color, linewidth=2)
        # Plot standard deviation as a shaded area
        plt.fill_between(x, mean_values - std_values, mean_values + std_values,
                         color=color, alpha=0.3, label=f"{os.path.basename(csv_file)} Std")

    file_names = []
    subject_ids = []
    for file in csv_files:
        file_names.append(os.path.basename(file))
    for dir in dirs:
        subject_ids.append(os.path.basename(dir))

    plt.xlabel("Along-Tract Node Position")
    plt.ylabel("FA Value")
    plt.title(f"FA Profiles for {file_names[0]}, subject {subject_ids[0]} "
              f"and {file_names[1]}, subject {subject_ids[1]}")
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.show()


