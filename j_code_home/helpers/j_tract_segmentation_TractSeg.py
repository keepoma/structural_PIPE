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
    peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")

    output_dir = paths["tractseg_dir"]
    os.makedirs(output_dir, exist_ok=True)

    bundleseg_dir =  os.path.join(output_dir, "bundle_segmentations")
    endingsseg_dir =  os.path.join(output_dir, "endings_segmentations")
    tom_dir =  os.path.join(output_dir, "TOM")

    # Bundle segmentation
    if (not os.path.exists(bundleseg_dir)) or \
            (len([f for f in os.listdir(bundleseg_dir)
                  if os.path.isfile(os.path.join(bundleseg_dir, f))]) != 72):
        run_cmd([
            "TractSeg",
            "-i", peaks_path_individual_RF,
            "-o", output_dir,
            "--output_type", "tract_segmentation"
        ])
    else:
        print(f"{bundleseg_dir} has 72 files, skipping")

    # Endings segmentation
    if (not os.path.exists(endingsseg_dir)) or \
            (len([f for f in os.listdir(endingsseg_dir)
                  if os.path.isfile(os.path.join(endingsseg_dir, f))]) != 144):
        run_cmd([
            "TractSeg",
            "-i", peaks_path_individual_RF,
            "-o", output_dir,
            "--output_type", "endings_segmentation"
        ])
    else:
        print(f"{endingsseg_dir} has 144 files, skipping")

    # Tract Orientation Maps (TOM)
    if (not os.path.exists(tom_dir)) or \
            (len([f for f in os.listdir(tom_dir)
                  if os.path.isfile(os.path.join(tom_dir, f))]) != 72):
        run_cmd([
            "TractSeg",
            "-i", peaks_path_individual_RF,
            "-o", output_dir,
            "--output_type", "TOM"
        ])
    else:
        print(f"{tom_dir} has 72 files, skipping")


def tractography_resample_and_extract_metrics(paths, nthreads):
    """
    Process a single subject by looping through each tract.
    """

    # Build the full path to the tract_name.txt file
    tract_names_file = "/media/nas/nikita/j_code_f/helpers/tract_name.txt"
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
        wmrf_norm = os.path.join(paths["dwi_dir"], "wm_norm.mif")

        files = [f for f in os.listdir(tracking_dir) if os.path.isfile(os.path.join(tracking_dir, f))]

        # Track generation using tckgen
        if len(files) != 144:
            run_cmd([
                "tckgen",
                "-algorithm", "iFOD2",
                wmrf_norm, tck_path,
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
        else:
            print(f"skipping iFOD2 and resampling generation")

        # Sample ADC values along original and resampled tract
        if not os.path.exists(adc_csv):
            adc_mif = os.path.join(paths["dwi_dir"], "adc.mif")
            run_cmd([
                "tcksample",
                tck_path,
                adc_mif,
                adc_csv,
                "-force"
            ])


        if not os.path.exists(adc_n100_csv):
            run_cmd([
                "tcksample",
                tck_n100_path,
                adc_mif,
                adc_n100_csv,
                "-force"
            ])


        # Sample FA values along the tract
        fa_mif = os.path.join(paths["dwi_dir"], "fa.mif")

        if not os.path.exists(fa_csv):
            run_cmd([
                "tcksample",
                tck_path,
                fa_mif,
                fa_csv,
                "-force"
            ])


        if not os.path.exists(fa_n100_csv):
            run_cmd([
                "tcksample",
                tck_n100_path,
                fa_mif,
                fa_n100_csv,
                "-force"
            ])


        # Peaks along tract
        peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")

        if not os.path.exists(peaks_txt):
            run_cmd([
                "tcksample",
                tck_path,
                peaks_path_individual_RF,
                peaks_txt,
                "-force"
            ])


        if not os.path.exists(peaks_n100_txt):
            run_cmd([
                "tcksample",
                tck_n100_path,
                peaks_path_individual_RF,
                peaks_n100_txt,
                "-force"
            ])



def tractseg_tracking_and_tractometry(root, paths, subj_dir, session_dir, session):
    """
    Implementing tractsegts own Tractometry as an alternative/conjunction to iFOD2
    """

    peaks_path_individual_RF = os.path.join(paths["dwi_dir"], "fod_peaks_individual_RF.nii.gz")
    tractseg_output = os.path.join(paths["tractseg_dir"])
    tom_trackings = os.path.join(tractseg_output, "TOM_trackings")
    if not os.path.exists(tom_trackings):
        run_cmd([
            "Tracking",
            "-i", peaks_path_individual_RF,
            "-o", tractseg_output,
            "--nr_fibers", "5000"
        ])
    else:
        print(f"{tom_trackings} exists, assuming completing and skipping")

    subject_id = os.path.basename(subj_dir)
    tractseg_tractometry_dir = os.path.join(root, "group_analysis", session, "tractseg_tractometry")
    os.makedirs(tractseg_tractometry_dir, exist_ok=True)
    subj_tractseg_tractometry_dir = os.path.join(tractseg_tractometry_dir, subject_id)
    os.makedirs(subj_tractseg_tractometry_dir, exist_ok=True)
    tractseg_tractometry_output = os.path.join(subj_tractseg_tractometry_dir, f"Tractometry.csv")
    bundle_e_path = os.path.join(session_dir, "tractseg_output", "endings_segmentations")
    fa_nii = os.path.join(paths["dwi_dir"], "fa.nii.gz")
    tom_trackings = os.path.join(tractseg_output, "TOM_trackings")
    files = [f for f in os.listdir(subj_tractseg_tractometry_dir) if os.path.isfile(os.path.join(subj_tractseg_tractometry_dir, f))]
    if not os.path.exists(tractseg_tractometry_output) or len(files) != 1:
        run_cmd([
            "Tractometry",
            "-i", tom_trackings,
            "-o", tractseg_tractometry_output,
            "-e", bundle_e_path,
            "-s", fa_nii
        ])
    else:
        print(f"{tractseg_tractometry_output} exists, skipping")
