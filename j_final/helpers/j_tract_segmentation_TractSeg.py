import os
import pandas as pd
from .j_helpers import run_cmd

"""
This code is also able to run as a standalone if the preprocessing part is
to be skipped
"""


def tract_and_endings_segmentation_TOMs(paths):
    """
    Runs tract segmentation, endings segmentation, and generates tract orientation maps.
    """

    # Helper lambda for paths in the 5_dwi folder
    peaks_path_group_rf = os.path.join(paths["dwi_dir"], "fod_peaks_group_RF.nii.gz")

    output_dir = paths["tractseg_dir"]
    os.makedirs(output_dir, exist_ok=True)

    # Tract segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path_group_rf,
        "-o", output_dir,
        "--output_type", "tract_segmentation"
    ])

    # Endings segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path_group_rf,
        "-o", output_dir,
        "--output_type", "endings_segmentation"
    ])

    #  Tract Orientation Maps
    run_cmd([
        "TractSeg",
        "-i", peaks_path_group_rf,
        "-o", output_dir,
        "--output_type", "TOM"
    ])


def tractography_resample_and_extract_metrics(paths, nthreads):
    """
    Process a single subject by looping through each tract.
    """

    # Build the full path to the tract_name.txt file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tract_names_file = os.path.join(script_dir, "helpers", "tract_name.txt")
    tract_names = pd.read_csv(tract_names_file, header=None)[0].tolist()

    output_dir = paths["tractseg_dir"]

    for tract_name in tract_names:

        # Build file paths for segmentation files and outputs
        bundle_path = os.path.join(output_dir, "bundle_segmentations",
                                   f"{tract_name}.nii.gz")
        bundle_b_path = os.path.join(output_dir, "endings_segmentations",
                                     f"{tract_name}_b.nii.gz")
        bundle_e_path = os.path.join(output_dir, "endings_segmentations",
                                     f"{tract_name}_e.nii.gz")

        # Output tracking paths
        tracking_dir = os.path.join(output_dir, "FOD_iFOD2_trackings")
        os.makedirs(tracking_dir, exist_ok=True)
        tck_path = os.path.join(tracking_dir, f"{tract_name}.tck")
        tck_n100_path = os.path.join(tracking_dir, f"{tract_name}_N100.tck")

        # Along-tract output directory with 3 subdirectories for each metric
        along_dir = paths["at_dir"]
        os.makedirs(along_dir, exist_ok=True)
        adc_dir = os.path.join(paths["at_dir"], "ADC")
        os.makedirs(adc_dir, exist_ok=True)
        fa_dir = os.path.join(paths["at_dir"], "FA")
        os.makedirs(fa_dir, exist_ok=True)
        peaks_dir = os.path.join(paths["at_dir"], "peaks")
        os.makedirs(peaks_dir, exist_ok=True)
        adc_csv = os.path.join(adc_dir, f"{tract_name}_adc.csv")
        adc_n100_csv = os.path.join(adc_dir, f"{tract_name}_n100_adc.csv")
        fa_csv = os.path.join(fa_dir, f"{tract_name}_fa.csv")
        fa_n100_csv = os.path.join(fa_dir, f"{tract_name}_n100_fa.csv")
        peaks_txt = os.path.join(peaks_dir, f"{tract_name}_peaks.txt")
        peaks_n100_txt = os.path.join(peaks_dir, f"{tract_name}_n100_peaks.txt")
        wmrf_group_norm = os.path.join(paths["dwi_dir"], "wm_group_average_based_norm.mif")

        # Track generation using tckgen
        run_cmd([
            "tckgen",
            "-algorithm", "iFOD2",
            wmrf_group_norm, tck_path,
            "-seed_image", bundle_path,
            "-seed_unidirectional",
            "-mask", bundle_path,
            "-include", bundle_b_path,
            "-include", bundle_e_path,
            "-minlength", "40",
            "-maxlength", "250",
            "-seeds", "1000000",
            "-select", "2000",
            "-cutoff", "0.1",
            "-nthreads", str(nthreads),
            "-force"
        ])

        # Resampling using tckresample
        run_cmd([
            "tckresample",
            tck_path,
            "-num_points", "100",
            tck_n100_path,
            "-nthreads", str(nthreads),
            "-force"
        ])

        # Sample ADC values along original and resampled tract
        adc_mif = os.path.join(paths["dwi_dir"], "adc.mif")
        run_cmd([
            "tcksample",
            tck_path,
            adc_mif,
            adc_csv,
            "-force"
        ])

        run_cmd([
            "tcksample",
            tck_n100_path,
            adc_mif,
            adc_n100_csv,
            "-force"
        ])

        # Sample FA values along the tract
        fa_mif = os.path.join(paths["dwi_dir"], "fa.mif"),

        run_cmd([
            "tcksample",
            tck_path,
            fa_mif,
            fa_csv,
            "-force"
        ])

        run_cmd([
            "tcksample",
            tck_n100_path,
            fa_mif,
            fa_n100_csv,
            "-force"
        ])

        # Peaks along tract
        peaks_path_group_RF = os.path.join(paths["dwi_dir"], "fod_peaks_group_RF.nii.gz")

        run_cmd([
            "tcksample",
            tck_path,
            peaks_path_group_RF,
            peaks_txt,
            "-force"
        ])

        run_cmd([
            "tcksample",
            tck_n100_path,
            peaks_path_group_RF,
            peaks_n100_txt,
            "-force"
        ])


def tractseg_tracking_and_tractometry(root, paths, session_dir):
    """
    Implementing tractsegts own Tractometry as an alternative/conjunction to iFOD2
    """

    peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")
    peaks_path_group_RF = os.path.join(paths["dwi_dir"], "fod_peaks_group_RF.nii.gz")
    tractseg_output = os.path.join(paths["tractseg_dir"])

    run_cmd([
        "Tracking",
        "-i", peaks_path_group_RF,
        "-o", tractseg_output,
        "--nr_fibers", "5000"
    ])

    subject_id = os.path.basename(session_dir)
    tractseg_tractometry_dir = os.path.join(root, "group_analysis", "tractseg_tractometry")
    os.makedirs(tractseg_tractometry_dir, exist_ok=True)
    subj_tractseg_tractometry_dir = os.path.join(tractseg_tractometry_dir, subject_id)
    os.makedirs(subj_tractseg_tractometry_dir, exist_ok=True)
    tom_trackings = os.path.join(tractseg_output, "TOM_trackings")
    tractseg_tractometry_output = os.path.join(subj_tractseg_tractometry_dir, f"Tractometry.csv")
    bundle_e_path = os.path.join(session_dir, "tractseg_output", "endings_segmentations")
    fa_output_nii = os.path.join(paths["dwi_dir"], "fa.nii.gz")

    run_cmd([
        "Tractometry",
        "-i", tom_trackings,
        "-o", tractseg_tractometry_output,
        "-e", bundle_e_path,
        "-s", fa_output_nii
    ])
