import os
import preprocessing_functions as preproc
import statistical_analysis as sa
from helpers.helpers import (get_subject_paths, get_subject_dirs, ask_yes_no,
                             fancy_print, prompt_for_folder)
from registration import register_t1_and_5tt_to_dwi


"""
This script contains the function with the main pipeline that shares the code necessary
for both pipes. It's supposed to be imported and built upon.
"""


def preprocessing_pipeline(root, nthreads, do_hsvs):
    subject_dirs = get_subject_dirs(root)

    is_preprocessed = ask_yes_no("Is every subject in this folder preprocessed?")
    has_registration = ask_yes_no("Has the registration of T1 and 5tt to dwi been done?")

    if not is_preprocessed:
        # Pending to save these to config file or something in case of script re-run
        print("Please provide the following folder names: ")
        t1_folder = prompt_for_folder("006_T1w_MPR", "T1 scan")
        t2_folder = prompt_for_folder("008_T2w_SPC", "T2 scan")
        t2_df_folder = prompt_for_folder("009_t2_space_dark-fluid_sag_p2_iso_0_8", "T2 FLAIR")
        dwi_ap_folder = prompt_for_folder("016_dMRI_dir98_AP", "dMRI AP scan")
        dwi_pa_folder = prompt_for_folder("019_dMRI_dir98_PA", "dMRI PA scan")

    for subj_dir in subject_dirs:
        # Retrieve standard paths
        paths = get_subject_paths(subj_dir)
        fancy_print("Executing script for Subject:", subj_dir)

        if not is_preprocessed:
            fancy_print("Preprocessing", subj_dir)
            fancy_print("Converting Scans", subj_dir)
            preproc.convert_scans(paths, nthreads, t1_folder, t2_folder, t2_df_folder, dwi_ap_folder, dwi_pa_folder)
            fancy_print("Preprocessing dMRI Data", subj_dir)
            preproc.preprocess_dwi(paths, nthreads)
            fancy_print("Calculating Response Function", subj_dir)
            preproc.response_function(paths, nthreads)

    if len(subject_dirs) > 1:
        group_output_directory = os.path.join(root, "group_analysis")
        print(f"\n========= Calculating Group Response Function =========\n")
        sa.compute_group_response_functions(root, group_output_directory, nthreads)

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        if not is_preprocessed:
            fancy_print("Performing FOD and normalization", subj_dir)
            preproc.FOD_normalization_peaks(paths, root, nthreads)
    print(f"\n========= Preprocessing of all subjects COMPLETE =========\n")

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        if not has_registration:
            fancy_print("Registering T1 and 5tt to dMRI Space", subj_dir)
            register_t1_and_5tt_to_dwi(paths, nthreads, do_hsvs)
    print(f"\n========= Transformation matrix and registration of T1 and 5tt to DWI for all subjects COMPLETE =========\n")
