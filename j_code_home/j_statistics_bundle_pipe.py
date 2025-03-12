import os
import logging
import numpy as np
from helpers.j_helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from stat_helpers.j_regional_analysis import (run_fa_processing, process_arcuate_mean_fa,
                                              groupwise_af_ttest)



"""
Bundle Statistics pipe
"""


def main():
    args = get_args()
    subject_dirs = get_subject_dirs(args.root)
    logs(args.root)
    sessions = ['ses_pre', 'ses_post']
    for subj_dir in subject_dirs:
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)

            process_arcuate_mean_fa(paths, subj_dir)
            run_fa_processing(paths)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")
    groupwise_af_ttest(args.root, subject_dirs)

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise
