import glob
import os 
import pandas as pd

root = 'subjects'
tract_names = pd.read_csv('tract_name.txt', header=None)[0]

for subject_path in glob.glob(root+'/*'):
	for tract_name in tract_names:
		bundle_path = subject_path + '/tractseg_output/bundle_segmentations/' + tract_name + '.nii.gz'
		bundle_b_path = subject_path + '/tractseg_output/endings_segmentations/' + tract_name + '_b.nii.gz'
		bundle_e_path = subject_path + '/tractseg_output/endings_segmentations/' + tract_name + '_e.nii.gz'

		os.system('tckgen -algorithm iFOD2 '+subject_path+ '/5_DWI/wm.mif '+subject_path+'/tractseg_output/FOD_iFOD2_trackings/'+tract_name+'.tck -seed_image  '+bundle_path+' -mask '+bundle_path+' -include '+bundle_b_path+' -include '+bundle_e_path+ ' -minlength 40 -maxlength 250 -seeds 1000000 -select 2000 -cutoff 0.05')
		os.system('tckresample '+subject_path+'/tractseg_output/FOD_iFOD2_trackings/'+tract_name+'.tck -num_points 100 '+subject_path+'/tractseg_output/FOD_iFOD2_trackings/'+tract_name+'_N100.tck')
		os.system('tcksample '+subject_path+'/tractseg_output/FOD_iFOD2_trackings/'+tract_name+'_N120.tck '+ subject_path+'/5_DWI/adc.mif '+subject_path+'/along_tract/'+tract_name+'_adc.csv -fo')
		os.system('tcksample '+subject_path+'/tractseg_output/FOD_iFOD2_trackings/'+tract_name+'_N120.tck '+ subject_path+'/6_MRE/c_b0.nii '+subject_path+'/along_tract/'+tract_name+'_c.csv -fo')
		
