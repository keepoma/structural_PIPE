import os
import logging
from helpers.j_helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from helpers.j_atlas_generation import freesurfer_atlas_generation
from helpers.j_preprocessing_pipeline import preprocessing_pipeline
from helpers.j_registration import register_5tt_to_dwi
from helpers.j_tractography_and_TDIs import streamline_seeding, generate_tracks_and_sift
from helpers.j_connectome import connectome_generation, generate_weighted_connectome_matrices



"""
This is the main regional tractography pipeline. 
All the necessary functions are imported from other modules.
"""


def main():
    args = get_args()
    subject_dirs = get_subject_dirs(args.root)
    logs(args.root)
    preprocessing_pipeline(args.root, args.nthreads)
    sessions = ['ses_pre', 'ses_post']
    for subj_dir in subject_dirs:
        for session in sessions:
            try:
                session_dir = os.path.join(subj_dir, session)
                paths = get_subject_paths(session_dir)
                subj_ses = subj_dir + '_' + session
                fancy_print("Generating Freesurfer/HCP-based atlas", subj_ses)
                freesurfer_atlas_generation(paths, args.nthreads, subj_ses)
                fancy_print("registering 5tt to dwi", subj_ses)
                register_5tt_to_dwi(paths, args.nthreads)
                fancy_print("Performing streamline seeding", subj_ses)
                streamline_seeding(paths)
                fancy_print("Generating whole-brain tracks and applying SIFT", subj_ses)
                generate_tracks_and_sift(paths, args.nthreads)
                fancy_print("Generating connectome matrix", subj_ses)
                connectome_generation(paths, args.nthreads)
                fancy_print("Generating connectome weighting by metrics", subj_ses)
                generate_weighted_connectome_matrices(paths, args.nthreads)
            except Exception as error:
                logging.exception(f"An error occurred while processing {subj_dir} {session}: {error}")
                # Optionally continue to next iteration
                continue
        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise
