import os
import logging
from helpers.j_helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from helpers.j_preprocessing_pipeline import preprocessing_pipeline

from helpers.j_tract_segmentation_TractSeg import (tract_and_endings_segmentation_TOMs,
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
    preprocessing_pipeline(args.root, args.nthreads)
    sessions = ['ses_pre', 'ses_post']
    for subj_dir in subject_dirs:
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)
            subj_ses = subj_dir + '_' + session
            fancy_print("Segmenting tracts, endings, and generating TOMs", subj_ses)
            tract_and_endings_segmentation_TOMs(paths)
            fancy_print("iFOD2 Track generation, resampling and calculating metrics", subj_ses)
            tractography_resample_and_extract_metrics(paths, args.nthreads)
            fancy_print("TractSeg tracking and tractometry", subj_ses)
            tractseg_tracking_and_tractometry(args.root, paths, subj_dir, session_dir, session)


        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise
