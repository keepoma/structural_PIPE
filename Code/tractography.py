import glob
import os 
import pandas as pd

max_threads = str(max(1, os.cpu_count() - 10))

root = '/media/nas/nikita/sample_598'
tract_names = pd.read_csv('/home/nikita/Structural_PIPE/Code/tract_name.txt', header=None)[0]

# Iteration by subject
for subject_path in glob.glob(root+'/*'):
	# Iteration by tract
	for tract_name in tract_names:
		# Segmentation by bundle and start/end regions
		bundle_path = f"{subject_path}/2_nifti/tractseg_output/bundle_segmentations/{tract_name}.nii.gz"
		bundle_b_path = f"{subject_path}/2_nifti/tractseg_output/endings_segmentations/{tract_name}_b.nii.gz"
		bundle_e_path = f"{subject_path}/2_nifti/tractseg_output/endings_segmentations/{tract_name}_e.nii.gz"

		# Shell command execution
		# Track generation by bundle
		os.system(f"tckgen -algorithm iFOD2 {subject_path}/5_DWI/wm.mif "
				  f"{subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}.tck "
				  f"-seed_image {bundle_path} -mask {bundle_path} "
				  f"-include {bundle_b_path} -include {bundle_e_path} "
				  "-minlength 40 -maxlength 250 -seeds 1000000 -select 2000 -cutoff 0.05 "
				  f"-nthreads {max_threads}")

		# Resampling
		os.system(f"tckresample {subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}.tck "
				  f"-num_points 100 {subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}_N100.tck "
				  f"-nthreads {max_threads}")

		# Apparent Diffusion Coefficient
		os.system(f"tcksample {subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}_N100.tck "
				  f"{subject_path}/5_DWI/adc.mif "
				  f"{subject_path}/along_tract/{tract_name}_adc.csv -fo")

		# Fractional Anisotropy
		os.system(f"tcksample {subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}_N100.tck "
				  f"{subject_path}/5_DWI/fa.mif "
				  f"{subject_path}/along_tract/{tract_name}_fa.csv -fo")

		os.system(f"tcksample {subject_path}/2_nifti/tractseg_output/FOD_iFOD2_trackings/{tract_name}_N100.tck "
				  f"{subject_path}/2_nifti/b0.nii.gz "
				  f"{subject_path}/along_tract/{tract_name}_c.csv -fo")

# Review
# AF left TractSeg result
mrview "/media/nas/nikita/sample_598/raw/6_mif/t1_coreg.mif" -tractography.load "/media/nas/nikita/sample_598/raw/2_nifti/tractseg_output/FOD_iFOD2_trackings/AF_left_N100.tck" -mode 3
