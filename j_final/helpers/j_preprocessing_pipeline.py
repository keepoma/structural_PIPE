import os
from . import j_preprocessing_functions as preproc
from .j_helpers import (get_subject_paths, get_subject_dirs, fancy_print)
from .j_registration import register_t1_to_dwi


"""
This script contains the function with the preproc steps that are shared between both pipes.
Supposed to be imported and built upon.
"""


def preprocessing_pipeline(root, nthreads):
    subject_dirs = get_subject_dirs(root)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            subj_ses = subj_dir + ' ' + session
            fancy_print("Converting Scans", subj_ses)
            preproc.convert_scans(paths, nthreads)
            fancy_print("Preprocessing dMRI Data", subj_ses)
            preproc.preprocess_dwi(paths, nthreads)
            fancy_print("Calculating Tensors and dMRI metrics", subj_ses)
            preproc.calculate_tensors_and_dmri_metrics(paths, nthreads)
            fancy_print("Calculating Response Function", subj_ses)
            preproc.response_function(paths, nthreads)

    group_output_directory = os.path.join(root, "group_analysis")
    print(f"\n========= Calculating Group Response Function =========\n")
    preproc.compute_group_response_functions(root, group_output_directory, nthreads)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            subj_ses = subj_dir + ' ' + session_dir
            fancy_print("Performing FOD and normalization", subj_ses)
            preproc.FOD_normalization_peaks(paths, root, nthreads)

    for subj_dir in subject_dirs:
        for session in ['ses_pre', 'ses_post']:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            subj_ses = subj_dir +' '+ session_dir
            fancy_print("Registering T1 to dMRI Space", subj_ses)
            register_t1_to_dwi(paths, nthreads)