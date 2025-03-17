import os
import numpy as np
import pandas as pd
from helpers.j_helpers import (get_subject_paths_wo_ses, get_args,
                               get_subject_dirs, logs,
                               get_subject_paths)
from stat_helpers.j_regional_analysis import (run_metric_processing, process_track_mean,
                                              groupwise_track_ttest, groupwise_track_wilcoxon,
                                              plot_tractometry_results, groupwise_track_wilcoxon_non_resampled,
                                              get_nodewise_metric_range, report_nonresampled_metric_stats,
                                              convert_peaks_to_absolute)
from stat_helpers.j_regional_visualization import (plot_mean_metrics_along_nodes, plot_paired_mean_fa_non_resampled,
                                                   plot_histogram_mean_fa_non_resampled)



"""
Bundle Statistics pipe
"""


def main():
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    subject_dirs = get_subject_dirs(root)
    track_base = "CST"
    metrics = ["FA", "ADC", "peaks"]

    sessions = ['ses_pre']
    for subj_dir in subject_dirs:
        for session in sessions:
            session_dir = os.path.join(subj_dir, session)
            paths = get_subject_paths(session_dir)

            for metric in metrics:
                process_track_mean(paths, session_dir, track_base, metric)
                run_metric_processing(paths, session_dir, track_base, metric)
                if metric.lower() == "peaks":
                    # Convert peaks to absolute and generate _abs_statistics.csv
                    convert_peaks_to_absolute(paths, session_dir, track_base)

    for metric in metrics:

        groupwise_track_wilcoxon(root, subject_dirs, track_base, metric)

        print("non resampled stats: ")
        report_nonresampled_metric_stats(subject_dirs, track_base, metric)

    plot_mean_metrics_along_nodes(
        subject_dirs,
        track_base,
        metrics,
        alpha=0.7,
        save_figure=False
    )

    # Run group-level statistical tests for the resampled tracts
    #groupwise_track_ttest(root, subject_dirs, track_base, metric)

    """
    plot_paired_mean_fa_non_resampled(subject_dirs, track_base)
    plot_tractometry_results(root)
    plot_histogram_mean_fa_non_resampled(subject_dirs, track_base)
    """


if __name__ == '__main__':
    main()

