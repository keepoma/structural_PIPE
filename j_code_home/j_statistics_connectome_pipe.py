import os
import logging
import numpy as np
from helpers.j_helpers import (get_subject_paths, fancy_print,
                             get_args, get_subject_dirs,
                             logs)
from stat_helpers.j_connectome_analysis import (compute_connectivity_metrics,
                                                compute_metrics_for_weight_threshold_range, compute_metrics_for_prethresholded_matrix,
                                                find_top_nodes_by_strength, threshold_and_save_matrix_by_top_percent)



"""
Statistics pipe
"""


def process_connectivity_matrices(matrix, matrix_path, lookup_path, con_stats_dir):
    """
    Computes the 90th and 80th percentile thresholded matrices from the given matrix,
    then creates a list consisting of:
      - the unmodified matrix,
      - the 90th percentile thresholded matrix, and
      - the 80th percentile thresholded matrix.

    For each of these, it computes the connectivity metrics (using
    compute_metrics_for_prethresholded_matrix) and graphs the weight distribution.
    The results are saved in unique folders under the provided con_stats_dir.
    """

    # Obtain the thresholded matrices.
    matrix_t90 = threshold_and_save_matrix_by_top_percent(matrix, 90, matrix_path)
    matrix_t80 = threshold_and_save_matrix_by_top_percent(matrix, 80, matrix_path)

    # Create a list of tuples: (matrix, folder name, label suffix)
    matrices_list = [
        (matrix, "unthresholded", "unmodified"),
        (matrix_t90, "t90", "t90"),
        (matrix_t80, "t80", "t80")
    ]

    for (mat, folder_label, label_suffix) in matrices_list:
        # Create an output folder for this threshold level.
        folder_path = os.path.join(con_stats_dir, folder_label)
        os.makedirs(folder_path, exist_ok=True)

        # Define a base filename using the label suffix.
        base_filename = f"hcpmmp1_invleng_{label_suffix}"

        # Compute and save the connectivity metrics for the current matrix.
        print(f"Computing metrics for {label_suffix} matrix", base_filename)
        compute_metrics_for_prethresholded_matrix(
            mat,
            lookup_path,
            folder_path,
            base_filename,
            binarize=False,
            overwrite=True
        )


def main():
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    lookup_path = "/Users/nikitakaruzin/MRI/projects/BATMAN/Supplementary_Files/hcpmmp1_ordered.txt"
    subject_dirs = get_subject_dirs(root)

    for subj_dir in subject_dirs:
        sessions = ['ses_pre', 'ses_post']
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            subj_ses = subj_dir + '_' + session
            paths = get_subject_paths(session_dir)
            print(f"Working on {subj_ses}")


            # matrix csv path and matrix numpy generation
            matrix_path = os.path.join(paths["atlas_dir"], "hcpmmp1_scale_invlength.csv")
            matrix = np.genfromtxt(matrix_path, delimiter=',')

            output_dir = paths["con_stats_dir"]
            process_connectivity_matrices(matrix, matrix_path, lookup_path, output_dir)

            # compute_connectivity_metrics(hcpmmp1_minmax_t80, metrics_file_path)

        print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")

    print(f"\n========= PROCESSING COMPLETE FOR ALL SUBJECTS =========\n")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        logging.exception("An unexpected error occurred")
        raise

    """
       VISUALIZATION
    """
    # This is already in the paper, no need to do this for multiple subjects
    # visualize_matrix_weights(hcpmmp1_minmax_t80, ".")
    # visualize_matrix(hcpmmp1_minmax_t80, clip=False)
    # visualize_matrix_side_by_side(hcpmmp1_minmax_t80)
    # visualize_matrix_comparison
    # visualize_saved_metrics(threshold_to_node_csv, threshold_to_global_csv)
