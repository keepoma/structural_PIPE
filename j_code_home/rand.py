import subprocess
import glob
import os
import numpy as np

def extract_gradients(mif_file):
    """
    Run mrinfo with the -grad option on a .mif file and return an array of gradients.
    Each line is expected to contain four comma-separated values: x, y, z, b.
    """
    # Run mrinfo with -grad option
    cmd = ['mrinfo', '-grad', mif_file]
    try:
        output = subprocess.check_output(cmd, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {mif_file}: {e}")
        return None

    # Process the output: split into lines and then into values
    gradients = []
    for line in output.strip().splitlines():
        # Skip header lines if any (assuming each valid line contains a comma)
        if ',' in line:
            parts = line.split(',')
            if len(parts) != 4:
                continue
            try:
                # Convert to floats
                grad = list(map(float, parts))
                gradients.append(grad)
            except ValueError:
                continue
    return np.array(gradients)

def count_directions_per_shell(gradients, tol=50):
    """
    Count number of directions per shell using the 4th column (b-value)
    with a given tolerance for grouping.
    Volumes with b < tol are considered b0.
    """
    bvals = gradients[:, 3]
    shells = {}
    for b in bvals:
        if b < tol:
            key = 0.0
        else:
            # Check if this b-value is close to an existing shell
            matched = False
            for shell in shells.keys():
                if shell != 0.0 and abs(b - shell) < tol:
                    key = shell
                    matched = True
                    break
            if not matched:
                key = b
        shells[key] = shells.get(key, 0) + 1
    return shells

def main():
    # Adjust the base directory as needed
    base_dir = '/Users/nikitakaruzin/Desktop/Research/Picht/j_stats'
    # Assuming each subject is under sub-*, and each has ses_pre and ses_post folders.
    pattern = os.path.join(base_dir, 'sub-*', 'ses_*', 'dwi', 'dwi_den_unr_pre_unbia_newor.mif')
    files = glob.glob(pattern)
    print(f"Found {len(files)} files.")

    for mif_file in files:
        print(f"\nProcessing file: {mif_file}")
        gradients = extract_gradients(mif_file)
        if gradients is None or len(gradients) == 0:
            print("No gradient information found.")
            continue
        shells = count_directions_per_shell(gradients, tol=50)
        for shell, count in shells.items():
            if shell == 0.0:
                print(f"b≈0: {count} volumes")
            else:
                print(f"b≈{shell:.1f}: {count} directions")

if __name__ == '__main__':
    main()
