import os
import logging
from helpers.helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from preprocessing_pipeline import preprocessing_pipeline
from tract_segmentation_TractSeg import (tract_and_endings_segmentation_TOMs,
                                         tractography_resample_and_extract_metrics,
                                         tractseg_tracking_and_tractometry)


"""
This is the main regional tractography pipeline. 
All the necessary functions are imported from other modules.
"""


def main():
    args = get_args()
    subject_dirs = get_subject_dirs(args.root)
    logs(args.root)
    preprocessing_pipeline(args.root, args.nthreads, do_hsvs=False)

    for subj_dir in subject_dirs:
        paths = get_subject_paths(subj_dir)
        fancy_print("Segmenting tracts, endings, and generating TOMs", subj_dir)
        tract_and_endings_segmentation_TOMs(paths, subj_dir)
        fancy_print("iFOD2 Track generation, resampling and calculating metrics", subj_dir)
        tractography_resample_and_extract_metrics(subj_dir, args.nthreads)
        fancy_print("TractSeg tracking and tractometry", subj_dir)
        tractseg_tracking_and_tractometry(args.root, paths, subj_dir)
        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise