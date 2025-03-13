import os
import logging
import numpy as np
from helpers.j_helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from stat_helpers.j_connectome_analysis import (compute_connectivity_metrics,
                                                compute_metrics_for_weight_threshold_range,
                                                find_top_nodes_by_strength)


"""
Statistics pipe
"""


def main():
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    subject_dirs = get_subject_dirs(root)
    for subj_dir in subject_dirs:
        sessions = ['ses_pre', 'ses_post']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            subj_ses = subj_dir + '_' + session
            paths = get_subject_paths(session_dir)

            print(f"Working on {subj_ses}")
            file_path = os.path.join(paths["atlas_dir"], "hcpmmp1_minmax.csv")
            matrix = np.genfromtxt(file_path, delimiter=',')
            lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/Supplementary_Files/hcpmmp1_ordered.txt"

            fancy_print("Computing connectivity metrics", subj_ses)
            os.makedirs(paths["con_stats_dir"], exist_ok=True)
            metrics_file_path = os.path.join(paths["con_stats_dir"], "connectivity_metrics.csv")
            compute_connectivity_metrics(matrix, metrics_file_path)

            fancy_print("Computing metrics for threshold range", subj_ses)
            thresholds = np.linspace(0.0, 500, 5)
            compute_metrics_for_weight_threshold_range(
                paths=paths,
                sc_path=file_path,
                lookup_path=lookup_path,
                thresholds=thresholds,
                binarize=False,
                overwrite=False
            )

            fancy_print("Localizing top nodes by strength", subj_ses)
            nodes_file_path = os.path.join(paths["con_stats_dir"], "nodes_strength.csv")
            find_top_nodes_by_strength(matrix, lookup_path,
                                       None, nodes_file_path)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise
