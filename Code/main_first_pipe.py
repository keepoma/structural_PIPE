import os
import preprocess_MRI_data as preproc
import statistical_analysis as sa
from helpers.helpers import run_cmd, get_subject_paths, get_subject_dirs, get_args, ask_yes_no, fancy_print
from tractography_TractSeg import tractography_resample_and_extract_metrics
from registration import register_t1_and_5tt_to_dwi

"""
Main pipeline code. 
"""

def tractseg(paths, subject_dir):
    """
    Runs tract segmentation, endings segmentation, and generates tract orientation maps.
    """

    # Helper lambda for paths in the 5_dwi folder
    peaks_path = os.path.join(paths["two_nifti"], "fod_peaks_individual_RF.nii.gz")

    # Tract segmentation
    output_dir = os.path.join(subject_dir, "tractseg_output")
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
        "--output_type", "tract_segmentation"
    ])

    # Endings segmentation
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
        "--output_type", "endings_segmentation"
    ])

    #  Tract Orientation Maps
    run_cmd([
        "TractSeg",
        "-i", peaks_path,
        "-o", output_dir,
        "--output_type", "TOM"
    ])


def tensor_and_scalar_metrics(paths, nthreads):
    """
    Runs diffusion tensor and computes related metrics
    """

    five_path = lambda subpath: os.path.join(paths["five_dwi"], subpath)

    # Diffusion tensor and ADC/FA computation
    run_cmd([
        "dwi2tensor", "-nthreads", str(nthreads),
        five_path("dwi_den_unr_pre_unbia.mif"),
        five_path("tensor.mif"), "-force"
    ])

    run_cmd([
        "tensor2metric", "-nthreads", str(nthreads),
        five_path("tensor.mif"),
        "-adc", five_path("adc.mif"), "-force"
    ])

    run_cmd([
        "tensor2metric", "-nthreads", str(nthreads),
        five_path("tensor.mif"),
        "-fa", five_path("fa.mif"), "-force"
    ])


def main():
    args = get_args()
    root = os.path.abspath(args.root)
    subject_dirs = get_subject_dirs(root, "group_analysis")

    is_preprocessed = ask_yes_no("Is every subject in this folder preprocessed?")
    has_registration = ask_yes_no("Has the registration of T1 and 5tt to dwi been done?")

    for subj_dir in subject_dirs:
        # Retrieve standard paths
        paths = get_subject_paths(subj_dir)
        fancy_print("Executing script for Subject:", subj_dir)

        if not is_preprocessed:
            fancy_print("Preprocessing", subj_dir)
            fancy_print("Converting Scans", subj_dir)
            preproc.convert_scans(paths, args.nthreads)
            fancy_print("Preprocessing dMRI Data", subj_dir)
            preproc.preprocess_dwi(paths, args.nthreads)
            fancy_print("Calculating Response Function", subj_dir)
            preproc.response_function(paths, args.nthreads)

    if len(subject_dirs) > 1:
        group_output_directory = os.path.join(root, "group_analysis")
        print(f"\n========= Calculating Group Response Function =========\n")
        sa.compute_group_response_functions(root, group_output_directory, args.nthreads)

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)

        if not is_preprocessed:
            fancy_print("Performing FOD and normalization", subj_dir)
            preproc.FOD_normalization_peaks(paths, root, args.nthreads)

        fancy_print("Running tractography", subj_dir)
        tractseg(paths, subj_dir)

        fancy_print("Generating Tensor and Scalar Metrics", subj_dir)
        tensor_and_scalar_metrics(paths, args.nthreads)

        if not has_registration:
            fancy_print("Registering T1 to dMRI Space", subj_dir)
            register_t1_and_5tt_to_dwi(paths, args.nthreads)

        fancy_print("Track generation, resampling and metrics generation", subj_dir)
        tractography_resample_and_extract_metrics(subj_dir)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")



if __name__ == "__main__":
    main()
