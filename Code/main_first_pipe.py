import os
import logging
from helpers.helpers import (get_subject_paths, fancy_print,
                             calculate_tensors_and_dmri_metrics,
                             get_args, get_subject_dirs, logs)
from tract_segmentation_TractSeg import (tract_and_endings_segmentation_TOMs,
                                         tractography_resample_and_extract_metrics)
from preprocessing_pipeline import general_pipeline


"""
This is the main regional tractography pipeline. 
All the necessary functions are imported from other modules.
"""


def main():
    args = get_args()
    subject_dirs = get_subject_dirs(args.root)
    logs(args.root)
    general_pipeline(args.root, args.nthreads)
    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        fancy_print("Running tractography", subj_dir)
        tract_and_endings_segmentation_TOMs(paths, subj_dir)
        fancy_print("Generating Tensor and Scalar Metrics", subj_dir)
        calculate_tensors_and_dmri_metrics(paths, args.nthreads)
        fancy_print("Track generation, resampling and metrics generation", subj_dir)
        tractography_resample_and_extract_metrics(subj_dir)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise