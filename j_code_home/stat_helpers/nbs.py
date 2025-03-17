import os
import subprocess
import numpy as np

def run_cmd(cmd, cwd=None):
    """Runs a system command and prints it."""
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)

def get_subject_dirs(root):
    """
    Returns a sorted list of subject directories in the given root.
    Assumes subject directories start with 'sub-'.
    """
    return sorted([
        os.path.join(root, d)
        for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)) and d.startswith("sub-")
    ])

def main():
    root = "/Users/nikitakaruzin/Desktop/Research/Picht/j_stats"
    subject_dirs = get_subject_dirs(root)

    pre_paths = []
    post_paths = []

    group_dir = os.path.join(root, "group_analysis")
    os.makedirs(group_dir, exist_ok=True)
    """


    # Loop through subjects and collect pre and post CSV files
    for subj_dir in subject_dirs:
        pre_file = os.path.join(subj_dir, "ses_pre", "atlas", "hcpmmp1_minmax.csv")
        post_file = os.path.join(subj_dir, "ses_post", "atlas", "hcpmmp1_minmax.csv")
        pre_paths.append(pre_file)
        post_paths.append(post_file)



    # 1) Create the input file for connectomes (connectomes.txt)
    connectomes_file = os.path.join(group_dir, "connectomes.txt")
    with open(connectomes_file, "w") as f:
        # Pre first
        for path in pre_paths:
            rel_path = os.path.relpath(path, start=group_dir)
            f.write(rel_path + "\n")
        # Then post
        for path in post_paths:
            rel_path = os.path.relpath(path, start=group_dir)
            f.write(rel_path + "\n")
    print(f"Connectomes file written to: {connectomes_file}")

    # 2) Create a 2-column design matrix:
    #    Column 1 = intercept (always 1)
    #    Column 2 = condition (0 for pre, 1 for post)
    design_file = os.path.join(group_dir, "designmatrix.txt")
    with open(design_file, "w") as f:
        # For all pre sessions
        for _ in pre_paths:
            f.write("1 0\n")  # [intercept=1, condition=0]
        # For all post sessions
        for _ in post_paths:
            f.write("1 1\n")  # [intercept=1, condition=1]
    print(f"Design matrix file written to: {design_file}")

    # 3) Create a contrast file with 2 entries in one line: [0 1]
    #    This ignores the intercept and tests the condition column
    contrast_file = os.path.join(group_dir, "contrast.txt")
    with open(contrast_file, "w") as f:
        f.write("0 1\n")
    print(f"Contrast file written to: {contrast_file}")

    # 4) Define an output prefix for connectomestats
    output_prefix = os.path.join(group_dir, "tfnbs_output")

    # 5) Run connectomestats
    cmd = [
        "connectomestats",
        connectomes_file,
        "tfnbs",   # algorithm
        design_file,
        contrast_file,
        output_prefix
    ]
    run_cmd(cmd, cwd=group_dir)
    print("connectomestats finished.")

    
    """

    # ---------------------------------------------------------
    # OPTIONAL: Analyze the final output (if you want to load
    # the corrected p-values & effect sizes).
    # Adjust the file names if they differ, e.g. tfnbs_outputfwe_1mpvalue_t2.csv
    # and tfnbs_outputbeta_1.csv.
    # ---------------------------------------------------------

    pvals_1mp_file = os.path.join(group_dir, "tfnbs_outputfwe_1mpvalue.csv")
    effects_file   = os.path.join(group_dir, "tfnbs_outputbeta_0.csv")

    if os.path.exists(pvals_1mp_file) and os.path.exists(effects_file):
        # Load (1 - p), convert to p
        print(f"Loading 1 - p-values from: {pvals_1mp_file}")
        vals = np.loadtxt(pvals_1mp_file, delimiter=",")
        pvals = 1.0 - vals

        print(f"Loading effect sizes (beta_1) from: {effects_file}")
        effects = np.loadtxt(effects_file, delimiter=",")

        if pvals.shape != effects.shape:
            raise ValueError(
                f"Shape mismatch!\n"
                f"pvals: {pvals.shape}\n"
                f"effects: {effects.shape}"
            )

        alpha = 0.05
        significant_mask = (pvals < alpha)
        binary_adjacency = significant_mask.astype(int)
        sig_effects = np.where(significant_mask, effects, 0)

        total_connections = pvals.size
        num_sig = significant_mask.sum()
        print(f"\n--- TFNBS Analysis Summary ---")
        print(f"Total connections: {total_connections}")
        print(f"Significant connections (p < {alpha}): {num_sig}")

        if num_sig > 0:
            sig_vals = effects[significant_mask]
            min_eff, max_eff = sig_vals.min(), sig_vals.max()
            print(f"Effect size range (significant only): {min_eff:.4f} to {max_eff:.4f}")
        else:
            print("No connections were significant at the chosen alpha level.")

        adjacency_out = os.path.join(group_dir, "significant_adjacency.csv")
        np.savetxt(adjacency_out, binary_adjacency, delimiter=",", fmt="%d")

        effects_out = os.path.join(group_dir, "significant_effects.csv")
        np.savetxt(effects_out, sig_effects, delimiter=",", fmt="%.6f")

        print(f"\nBinary adjacency matrix saved to: {adjacency_out}")
        print(f"Masked effect size matrix saved to: {effects_out}")
    else:
        print("Did not find the expected TFNBS output files for post-analysis.")

if __name__ == '__main__':
    main()
