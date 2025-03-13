import os
import numpy as np
import pandas as pd
from helpers.j_helpers import (get_subject_paths_wo_ses, get_args,
                               get_subject_dirs, logs)
from stat_helpers.j_regional_analysis import (run_fa_processing, process_track_mean_fa,
                                              groupwise_track_ttest, groupwise_track_wilcoxon,
                                              plot_tractometry_results, groupwise_track_wilcoxon_non_resampled,
                                              get_nodewise_fa_range, report_nonresampled_fa_stats)
from stat_helpers.j_regional_visualization import (plot_mean_fa_along_nodes, plot_paired_mean_fa_non_resampled,
                                                   plot_histogram_mean_fa_non_resampled)



"""
Bundle Statistics pipe
"""


def main():
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats/Bundles"
    subject_dirs = get_subject_dirs(root)
    # Define the track base (e.g., "AF" or "CST")
    track_base = "AF"

    for subj_dir in subject_dirs:
        paths = get_subject_paths_wo_ses(subj_dir)

        # Process non-resampled mean FA values
        process_track_mean_fa(paths, subj_dir, track_base)

        # Process resampled FA stats (per node)
        run_fa_processing(paths, subj_dir, track_base)

    plot_mean_fa_along_nodes(subject_dirs, track_base)

    # Run group-level statistical tests for the resampled tracts
    groupwise_track_ttest(root, subject_dirs, track_base)
    groupwise_track_wilcoxon(root, subject_dirs, track_base)

    groupwise_track_wilcoxon_non_resampled(root, subject_dirs, track_base)

    get_nodewise_fa_range(subject_dirs, track_base)
    print("non resampled stats: ")
    report_nonresampled_fa_stats(subject_dirs, track_base)
    """
    plot_paired_mean_fa_non_resampled(subject_dirs, track_base)
    plot_tractometry_results(root)
    plot_histogram_mean_fa_non_resampled(subject_dirs, track_base)
    """


if __name__ == '__main__':
    main()

