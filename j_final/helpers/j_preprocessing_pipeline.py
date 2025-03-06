import os
import j_preprocessing_functions as preproc
from j_helpers import (get_subject_paths, get_subject_dirs,
                             fancy_print, prompt_for_folder)
from j_registration import register_t1_and_5tt_to_dwi


"""
This script contains the function with the preproc steps that are shared between both pipes.
Supposed to be imported and built upon.
"""


def preprocessing_pipeline(root, nthreads, do_hsvs):
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