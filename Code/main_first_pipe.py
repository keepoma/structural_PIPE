import os
import preprocess_MRI_data as preproc
import statistical_analysis as sa
from helpers.helpers import (get_subject_paths, get_subject_dirs,
                             get_args, ask_yes_no,
                             fancy_print, tensor_and_scalar_metrics)
from tractography_TractSeg import tract_and_endings_segmentation_TOMs, tractography_resample_and_extract_metrics
from registration import register_t1_and_5tt_to_dwi

"""
This is the main pipeline code. 
All the necessary functions are imported from other modules.
"""


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
    tract_and_endings_segmentation_TOMs(paths, subj_dir)

    fancy_print("Generating Tensor and Scalar Metrics", subj_dir)
    tensor_and_scalar_metrics(paths, args.nthreads)

    if not has_registration:
        fancy_print("Registering T1 to dMRI Space", subj_dir)
        register_t1_and_5tt_to_dwi(paths, args.nthreads)

    fancy_print("Track generation, resampling and metrics generation", subj_dir)
    tractography_resample_and_extract_metrics(subj_dir)

    print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")
