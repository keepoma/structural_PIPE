# Structural Data Processing Pipelines (BIDS-Compatible)

This folder contains pipelines adapted to the BIDS format for structural data processing. Before running the pipelines, ensure all requirements are met.

## 1. Software Requirements

Make sure the following software packages are installed:
-MRtrix3
-ANTs
-FSL
-TractSeg

## 2. Python Version Compatibility
For proper functioning with MRtrix (which uses the outdated "imp" module), Python must be ≤ 3.11.
To ensure compatibility, create and activate a Conda environment with Python 3.11. 

## 3. Running the Pipelines
After activating your Conda environment execute the desired pipeline from the terminal:
-**j_complete_pipe.py**
Runs preprocessing and both the regional and connectome generation pipelines for pre- and post-therapy sessions for all subjects in the specified root directory.
-**j_connectome_pipe.py**
Runs preprocessing and generates the connectome.
-**j_gen_bundle.py**
Runs preprocessing and generates bundles.

Note: The statistical analysis pipelines have been adapted for my specific needs and may not be relevant for other users.

## 4. Usage Example
The command expects a `--root` argument that points to the study subjects directory. For example:
`python3 j_complete_pipe.py --root "/media/nas/nikita/test_study"`
