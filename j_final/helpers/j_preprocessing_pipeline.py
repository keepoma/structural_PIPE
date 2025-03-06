import os
import j_preprocessing_functions as preproc
from j_helpers import (get_subject_paths, get_subject_dirs, ask_yes_no,
                             fancy_print, prompt_for_folder)
from j_registration import register_t1_and_5tt_to_dwi


"""
This script contains the function with the preproc steps that are shared between both pipes.
Also contains a function for individual subject preprocessing
Supposed to be imported and built upon.
"""


def preprocessing_pipeline(root, nthreads, do_hsvs):
    """
    Main preprocessing pipeline.
    """

    subject_dirs = get_subject_dirs(root)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            subj_ses = subj_dir +' '+ session_dir
            fancy_print("Converting Scans", subj_ses)
            preproc.convert_scans(session_dir, nthreads)
            fancy_print("Preprocessing dMRI Data", subj_ses)
            preproc.preprocess_dwi(session_dir, nthreads)
            fancy_print("Calculating Tensors and dMRI metrics", subj_ses)
            preproc.calculate_tensors_and_dmri_metrics(session_dir, nthreads)
            fancy_print("Calculating Response Function", subj_ses)
            preproc.response_function(session_dir, nthreads)

    group_output_directory = os.path.join(root, "group_analysis")
    print(f"\n========= Calculating Group Response Function =========\n")
    preproc.compute_group_response_functions(root, group_output_directory, nthreads)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            subj_ses = subj_dir +' '+ session_dir
            fancy_print("Performing FOD and normalization", subj_ses)
            preproc.FOD_normalization_peaks(session_dir, root, nthreads)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            subj_ses = subj_dir +' '+ session_dir
            fancy_print("Registering T1 and 5tt to dMRI Space", subj_ses)
            register_t1_and_5tt_to_dwi(session_dir, nthreads, do_hsvs)


def individual_preproc(root, nthreads, do_hsvs):
    """
    Function for preprocessing on individual subject basis
    """

    subject_dirs = get_subject_dirs(root)
    # Print all available directories with their index
    print("Available subject directories:")
    for i, directory in enumerate(subject_dirs):
        print(f"{i}: {directory}")

    # Prompt the user to select one by index
    try:
        index = int(input("Enter the index number of the subject directory to use: "))
        if 0 <= index < len(subject_dirs):
            subj_dir = subject_dirs[index]
            print(f"Selected subject directory: {subj_dir}")
        else:
            print("Error: The index is out of range.")
    except ValueError:
        print("Error: Please enter a valid integer.")


    print("Please provide the following folder names: ")
    t1_folder = prompt_for_folder("006_T1w_MPR", "T1 scan")
    t2_folder = prompt_for_folder("008_T2w_SPC", "T2 scan")
    t2_df_folder = prompt_for_folder("009_t2_space_dark-fluid_sag_p2_iso_0_8", "T2 FLAIR")
    dwi_ap_folder = prompt_for_folder("016_dMRI_dir98_AP", "dMRI AP scan")
    dwi_pa_folder = prompt_for_folder("019_dMRI_dir98_PA", "dMRI PA scan")

    # Retrieve standard paths
    paths = get_subject_paths(subj_dir)
    fancy_print("Converting Scans", subj_dir)
    preproc.convert_scans(paths, nthreads, t1_folder, t2_folder, t2_df_folder, dwi_ap_folder, dwi_pa_folder)
    fancy_print("Preprocessing dMRI Data", subj_dir)
    preproc.preprocess_dwi(paths, nthreads)
    fancy_print("Calculating Tensors and dMRI metrics", subj_dir)
    preproc.calculate_tensors_and_dmri_metrics(paths, nthreads)
    fancy_print("Calculating Response Function", subj_dir)
    preproc.response_function(paths, nthreads)

    group_output_directory = os.path.join(root, "group_analysis")
    print(f"\n========= Calculating Group Response Function =========\n")
    preproc.compute_group_response_functions(root, group_output_directory, nthreads)

    fancy_print("Performing FOD and normalization", subj_dir)
    preproc.FOD_normalization_peaks(paths, root, nthreads)
    print(f"\n========= Preprocessing of all subjects COMPLETE =========\n")

    fancy_print("Registering T1 and 5tt to dMRI Space", subj_dir)
    register_t1_and_5tt_to_dwi(paths, nthreads, do_hsvs)
    print(
        f"\n========= Transformation matrix and registration of T1 and 5tt to DWI COMPLETE =========\n")


def main():
    individual_preproc("/media/nas/nikita/5p_test_controls", 60, False)


if __name__ == '__main__':
    main()