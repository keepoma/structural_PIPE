import os
import logging
from datetime import datetime
from helpers.helpers import (get_subject_paths, fancy_print,
                             calculate_tensors_and_dmri_metrics,
                             get_args, get_subject_dirs)
from atlas_generation import freesurfer_atlas_generation
from tractography_and_TDIs import streamline_seeding, generate_tracks_and_sift, generate_tdis
from connectome import connectome_generation, generate_weighted_connectome_matrices
from general_pipeline import general_pipeline


"""
This is the main pipeline leading to connectome generation.  
All the necessary functions are imported from other modules.
Example run: python3 main_second_pipe.py --root /home/nikita/Nikita_MRI
"""

def main():
    args = get_args()
    subject_dirs = get_subject_dirs(args.root)

    # Create the logs folder under args.root
    log_dir = os.path.join(args.root, "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Create a timestamped log filename
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = os.path.join(log_dir, f"log_{timestamp}.log")

    # Configure logging to write to the log file
    logging.basicConfig(
        filename=log_filename,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s'
    )

    general_pipeline(args.root, args.nthreads)
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        subject_id = os.path.basename(subj_dir)
        fancy_print("Performing streamline seeding", subj_dir)
        streamline_seeding(paths)
        fancy_print("Generating whole-brain tracks and applying SIFT", subj_dir)
        generate_tracks_and_sift(paths, args.nthreads)
        fancy_print("Generating TDIs and aligning T1", subj_dir)
        generate_tdis(paths, args.nthreads)
        fancy_print("Generating Freesurfer/HCP-based atlas", subj_dir)
        freesurfer_atlas_generation(paths, args.nthreads, subject_id)
        fancy_print("Generating connectome matrix", subj_dir)
        connectome_generation(paths, args.nthreads)
        fancy_print("Calculating Tensor and related metrics", subj_dir)
        calculate_tensors_and_dmri_metrics(paths, args.nthreads)
        fancy_print("Generating connectome weighting by metrics", subj_dir)
        generate_weighted_connectome_matrices(paths, args.nthreads)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise