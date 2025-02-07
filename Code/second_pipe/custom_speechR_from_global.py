import subprocess
import itertools


"""
Just looping through Inferior Frontal Gyrus and Superior Temporal Gyrus for now:
uses "44", "45", "47l", "STGa", "TA2", "STSdp", "STSda", "STSvp", "STSva"
generates 18 connections
"""

name_parcelindices = {
    "44": 74,
    "45": 75,
    "47l": 76,
    "IFJp": 77,
    "IFJa": 78,
    "IFSp": 79,
    "IFSa": 80,
    "STGa": 123,
    "TA2": 107,
    "STSdp": 129,
    "STSda": 128,
    "STSvp": 130,
    "STSva": 176,
    "TGd": 140,
    "TGv": 141
}

# Define base command components
tck_file = "/media/nas/nikita/sample_598/raw/8_tck/sift_1mio.tck"
assignments_file = "/media/nas/nikita/sample_598/raw/9_atlas/assignments_hcpmmp1.csv"

# Get a list of values from the dictionary
values = list(name_parcelindices.values())

# Loop over every unique pair of values (nodes)
for node_a, node_b in itertools.combinations(values, 2):
    # Create a custom output filename using the node values
    output_file = f"/media/nas/nikita/sample_598/raw/8_tck/speech_test/speech"
    command = [
        "connectome2tck",
        "-nodes", f"{node_a},{node_b}",
        "-exclusive",
        tck_file,
        assignments_file,
        output_file,
        "-nthreads", "102",
        "-force"
    ]

    print(f"Running: {' '.join(command)}")
    subprocess.run(command)

print("Processing complete!")
# then MERGED into 1 file using:
# tckedit /media/nas/nikita/sample_598/raw/8_tck/speech_test/speech*.tck /media/nas/nikita/sample_598/raw/8_tck/speech_test/left_speech_merged.tck -nthreads 102