import os
from helpers.helpers import (get_subject_paths, fancy_print,
                             calculate_tensors_and_dmri_metrics)
from tract_segmentation_TractSeg import (tract_and_endings_segmentation_TOMs,
                                         tractography_resample_and_extract_metrics)
from general_pipeline import general_pipeline


"""
This is the main regional tractography pipeline. 
All the necessary functions are imported from other modules.
"""


subject_dirs, args = general_pipeline()
for subj_dir in subject_dirs:
    paths = get_subject_paths(subj_dir)
    fancy_print("Running tractography", subj_dir)
    tract_and_endings_segmentation_TOMs(paths, subj_dir)
    fancy_print("Generating Tensor and Scalar Metrics", subj_dir)
    calculate_tensors_and_dmri_metrics(paths, args.nthreads)
    fancy_print("Track generation, resampling and metrics generation", subj_dir)
    tractography_resample_and_extract_metrics(subj_dir)

    print(f"\n========= Subject: {os.path.basename(subj_dir)} COMPLETE =========\n")
