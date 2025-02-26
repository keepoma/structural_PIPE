import os
import preprocess_MRI_data as preproc
import statistical_analysis as sa
from helpers.helpers import (run_cmd, get_subject_dirs,
                             get_subject_paths,get_args,
                             ask_yes_no, fancy_print,
                             calculate_tensors_and_dmri_metrics)
from registration import register_t1_and_5tt_to_dwi
from atlas_generation import freesurfer_atlas_generation
from tractography_and_TDIs import streamline_seeding, generate_tracks_and_sift, generate_tdis
from connectome import connectome_generation, generate_weighted_connectome_matrices

"""
Main pipeline for connectome construction
Example run: python3 main_second_pipe.py --root /home/nikita/Nikita_MRI

"""


args = get_args()
root = os.path.abspath(args.root)
subject_dirs = get_subject_dirs(root, "group_analysis")

is_preprocessed = ask_yes_no("Is every subject in this folder preprocessed?")
has_registration = ask_yes_no("Has the registration of T1 and 5tt to dwi been done?")

for subj_dir in subject_dirs:
    # Retrieve standard paths
    paths = get_subject_paths(subj_dir)
    subject_id = os.path.basename(subj_dir)
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

    if not has_registration:
        fancy_print("Registering T1 and 5tt to dMRI Space", subj_dir)
        register_t1_and_5tt_to_dwi(paths, args.nthreads)

    fancy_print("Performing streamline seeding", subj_dir)
    streamline_seeding(paths)
    fancy_print("Generating whole-brain tracks and applying SIFT", subj_dir)
    generate_tracks_and_sift(paths, args.nthreads)
    fancy_print("Generating TDIs and aligning T1", subj_dir)
    generate_tdis(paths, args.nthreads)
    # roi_localization(paths, args.nthreads)
    fancy_print("Generating Freesurfer/HCP-based atlas", subj_dir)
    freesurfer_atlas_generation(paths, args.nthreads, subject_id)
    fancy_print("Generating connectome matrix", subj_dir)
    connectome_generation(paths, args.nthreads)
    fancy_print("Calculating Tensor and related metrics", subj_dir)
    calculate_tensors_and_dmri_metrics(paths, args.nthreads)
    fancy_print("Generating connectome weighting by metrics", subj_dir)
    generate_weighted_connectome_matrices(paths, args.nthreads)

    print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")