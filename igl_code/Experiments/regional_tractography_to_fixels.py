import os
from helpers import run_cmd, get_args, get_subject_dirs, get_subject_paths

def warp_and_combine_arcuate(root, nthreads, tck_in_subj, combined_tck):
    """
    1) For each subject, warp 'tck_in_subj' to template space using subject2template_warp.mif.
    2) Combine (concatenate) all warped .tck into a single 'combined_tck' in template space.
       This yields a group-level .tck file that can represent the union of all subjects' arcuate.
    """

    template_dir = os.path.join(root, "group_analysis", "template")
    combined_tck_path = os.path.join(template_dir, combined_tck)

    # First remove any existing combined file (if re-running).
    if os.path.exists(combined_tck_path):
        os.remove(combined_tck_path)

    subject_dirs = get_subject_dirs(root)
    warped_tck_list = []

    for subj_dir in subject_dirs:
        subject_id = os.path.basename(subj_dir)
        paths = get_subject_paths(subj_dir)
        warp_file = os.path.join(paths["five_dwi"], "subject2template_warp.mif")

        if not os.path.isfile(warp_file):
            print(f"[WARNING] No warp file for {subject_id} at {warp_file}. Skipping.")
            continue

        subject_tck_in = os.path.join(subj_dir, "tractseg_output", "FOD_iFOD2_trackings", tck_in_subj)
        if not os.path.isfile(subject_tck_in):
            print(f"[WARNING] No tract file for {subject_id} at {subject_tck_in}. Skipping.")
            continue

        # Warp .tck to template space
        warped_tck_name = f"{subject_id}_arcuate_in_template.tck"
        warped_tck_path = os.path.join(template_dir, warped_tck_name)
        cmd = [
            "tcktransform",
            subject_tck_in,
            warp_file,
            warped_tck_path,
            "-nthreads", str(nthreads),
            "-force"
        ]
        run_cmd(cmd)
        warped_tck_list.append(warped_tck_path)

    # Now combine all warped TCKs into one big file (union)
    # We'll just 'cat' them. Another approach is 'tckedit' with multiple inputs.
    if len(warped_tck_list) == 0:
        raise RuntimeError("No warped .tck files found. Check your input paths / subject data.")

    # tckedit can combine multiple TCKs:
    # e.g. tckedit subj1.tck subj2.tck ... combined_tck
    tckedit_cmd = [
        "tckedit",
        *warped_tck_list,
        combined_tck_path,
        "-force"
    ]
    run_cmd(tckedit_cmd)

    print(f"[INFO] Created combined group .tck: {combined_tck_path}\n"
          f"Included {len(warped_tck_list)} warped arcs.")


def run_tract_mask_fba(root, nthreads, group_tck_name, tract_prefix):
    """
    Using a SINGLE combined group .tck in template space, subset FD/FC/FDC to that tract,
    then smooth, then run fixelcfestats.
    """

    template_dir = os.path.join(root, "group_analysis", "template")
    fixel_mask_dir = os.path.join(template_dir, "fixel_mask")
    matrix_dir = os.path.join(template_dir, "matrix")
    wmfod_template = os.path.join(template_dir, "wmfod_template.mif")

    group_tck_path = os.path.join(template_dir, group_tck_name)
    if not os.path.isfile(group_tck_path):
        raise FileNotFoundError(f"Missing group TCK in template space: {group_tck_path}")

    # 1) tck2fixel -> get fixel-level mask for the entire group
    tract_fixel_mask_dir = os.path.join(template_dir, f"{tract_prefix}_fixel_mask")
    run_cmd([
        "tck2fixel",
        group_tck_path,
        fixel_mask_dir,
        tract_fixel_mask_dir,
        "mask.mif",
        "-force"
    ])

    # 2) Subset FD/log_fc/fdc
    fd_sub_dir   = os.path.join(template_dir, f"fd_{tract_prefix}")
    logfc_sub_dir= os.path.join(template_dir, f"log_fc_{tract_prefix}")
    fdc_sub_dir  = os.path.join(template_dir, f"fdc_{tract_prefix}")

    run_cmd([
        "fixelfilter",
        os.path.join(template_dir, "fd"),
        "select",
        fd_sub_dir,
        "-mask", tract_fixel_mask_dir,
        "-matrix", matrix_dir,
        "-force"
    ])
    run_cmd([
        "fixelfilter",
        os.path.join(template_dir, "log_fc"),
        "select",
        logfc_sub_dir,
        "-mask", tract_fixel_mask_dir,
        "-matrix", matrix_dir,
        "-force"
    ])
    run_cmd([
        "fixelfilter",
        os.path.join(template_dir, "fdc"),
        "select",
        fdc_sub_dir,
        "-mask", tract_fixel_mask_dir,
        "-matrix", matrix_dir,
        "-force"
    ])

    # 3) Smooth
    fd_sub_smooth_dir   = f"{fd_sub_dir}_smooth"
    logfc_sub_smooth_dir= f"{logfc_sub_dir}_smooth"
    fdc_sub_smooth_dir  = f"{fdc_sub_dir}_smooth"

    run_cmd([
        "fixelfilter",
        fd_sub_dir,
        "smooth",
        fd_sub_smooth_dir,
        "-matrix", matrix_dir,
        "-force"
    ])
    run_cmd([
        "fixelfilter",
        logfc_sub_dir,
        "smooth",
        logfc_sub_smooth_dir,
        "-matrix", matrix_dir,
        "-force"
    ])
    run_cmd([
        "fixelfilter",
        fdc_sub_dir,
        "smooth",
        fdc_sub_smooth_dir,
        "-matrix", matrix_dir,
        "-force"
    ])

    # 4) fixelcfestats
    files_list = os.path.join(template_dir, "files.txt")
    design_txt = os.path.join(template_dir, "design_matrix.txt")
    contrast_txt = os.path.join(template_dir, "contrast_matrix.txt")

    stats_fd_dir   = os.path.join(template_dir, f"stats_fd_{tract_prefix}")
    stats_logfc_dir= os.path.join(template_dir, f"stats_log_fc_{tract_prefix}")
    stats_fdc_dir  = os.path.join(template_dir, f"stats_fdc_{tract_prefix}")
    os.makedirs(stats_fd_dir, exist_ok=True)
    os.makedirs(stats_logfc_dir, exist_ok=True)
    os.makedirs(stats_fdc_dir, exist_ok=True)

    run_cmd([
        "fixelcfestats",
        fd_sub_smooth_dir,
        files_list,
        design_txt,
        contrast_txt,
        matrix_dir,
        stats_fd_dir,
        "-force"
    ])
    run_cmd([
        "fixelcfestats",
        logfc_sub_smooth_dir,
        files_list,
        design_txt,
        contrast_txt,
        matrix_dir,
        stats_logfc_dir,
        "-force"
    ])
    run_cmd([
        "fixelcfestats",
        fdc_sub_smooth_dir,
        files_list,
        design_txt,
        contrast_txt,
        matrix_dir,
        stats_fdc_dir,
        "-force"
    ])

    print(f"[INFO] Finished tract-specific FBA on group mask: {tract_prefix}\n"
          f"  Fixel mask: {tract_fixel_mask_dir}\n"
          f"  FD subset:  {fd_sub_dir}\n"
          f"  log_fc subset: {logfc_sub_dir}\n"
          f"  fdc subset: {fdc_sub_dir}\n"
          f"  Stats FD:   {stats_fd_dir}\n"
          f"  Stats log FC: {stats_logfc_dir}\n"
          f"  Stats FDC:  {stats_fdc_dir}\n")


if __name__ == "__main__":
    args = get_args()
    root = os.path.abspath(args.root)

    # 1) Warp each subject's .tck to template space, then combine them into a single "arcuate_combined.tck"
    # This is optional if you already have a single .tck for the group in template space.
    warp_and_combine_arcuate(
        root=root,
        nthreads=args.nthreads,
        tck_in_subj="AF_left.tck",  # the .tck name each subject has in subject space
        combined_tck="AF_left_combined.tck"
    )

    # 2) Subset FD/FC/FDC using that combined group tract:
    run_tract_mask_fba(
        root=root,
        nthreads=args.nthreads,
        group_tck_name="AF_left_combined.tck",
        tract_prefix="AF_left"
    )
