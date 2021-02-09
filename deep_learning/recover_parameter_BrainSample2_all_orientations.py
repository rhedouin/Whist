
import os
import sys

import numpy as np
import scipy as sp

import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.io

import keras
import tensorflow

import h5py

import nibabel as nib


# Set parameters
nb_orientations = 9
nb_TE = 12
length_one_orientation = 2*nb_TE + 1

# Prepare Input name
data_project_folder = '/project/3015069.01/'
print(data_project_folder)
exp_name = 'BrainSample2LorentzinaCorrection'

list_fa = ['fa-05', 'fa-10', 'fa-15','fa-20', 'fa-35', 'fa-60']
nbFa =  len(list_fa)

noise = '4'

for kFa in range(0, nbFa):
	fa = list_fa[kFa]

	base_signal_name = 'BrainSample-2_ses-03'
	base_signal_folder = os.path.join(data_project_folder, 'derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/')
	mask_folder = os.path.join(data_project_folder, 'derived/BrainSample-2/ses-03/gre_renaud/masks/')
	print(base_signal_folder)

	dp_project_folder = '/project/3015069.04/'
	base_dp_folder = os.path.join(dp_project_folder, 'deep_learning/multi_orientations/BrainSample2LorentzinaCorrection/')
	print(base_dp_folder)

	dict_base_name = 'SignalWithNoise' + noise + '_8rep_' + str(nb_orientations) + 'rotations_' + str(nb_TE) + 'TE_' + exp_name + '_polyfit_cartesian_with_theta'
	dp_folder = os.path.join(base_dp_folder, dict_base_name + '/')
	model_path = os.path.join(dp_folder, 'epoch_40.h5py')

	input_scale_name = dict_base_name + '_input_scale.npy'
	input_scale_path = os.path.join(dp_folder, input_scale_name)

	input_signal_nii_path = os.path.join(base_signal_folder, base_signal_name + '_all_orientations_' + fa  +'_concatenate_signal_polyfit_cartesian_with_theta_2_ref.nii.gz')
	input_mask_nii_path = os.path.join(mask_folder, 'BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz')

	# Prepare Output name
	output_scale_min_path = os.path.join(dp_folder,dict_base_name + '_output_scale_min.npy')
	output_scale_max_path = os.path.join(dp_folder,dict_base_name + '_output_scale_max.npy') 
	
	output_suffix =    fa + '_polyfit_cartesian_with_theta_noise' +noise

	output_parameter_folder = os.path.join(base_signal_folder, 'parameter_maps/' + exp_name + '/' + fa)

	# Load/reshape data
	temp = nib.load(input_signal_nii_path)

	inputs = temp.get_fdata() 
	inputs = np.moveaxis(inputs, -1, 0)
	inputs_shape = inputs.shape

	image_length = inputs_shape[1] * inputs_shape[2] *inputs_shape[3]
	print("original image_length: ", image_length)

	inputs = inputs.reshape(inputs.shape[0],-1)
	inputs = np.transpose(inputs)

	# Load/reshape mask
	temp = nib.load(input_mask_nii_path)
	mask = temp.get_fdata()
	mask_shape = mask.shape
	print("mask shape: ", mask_shape)

	mask = mask.reshape(-1)
	mask = mask.astype(bool)

	# Set output nifti header
	affine = temp.affine
	header = temp.header
	header.set_data_dtype(16) 

	# Mask input
	inputs_masked = inputs[mask,:]
	print("masked input shape: ", inputs_masked.shape)  
	  
	# Load input/output normalization scale
	out_labels =  'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'

	input_scale = np.load(input_scale_path)
	output_scale_min = np.load(output_scale_min_path)
	output_scale_max = np.load(output_scale_max_path)

	len_outputs = len(output_scale_min)            
	len_inputs = inputs_shape[0]

	print("Len Inputs: ", len_inputs)
	print("Len Outputs: ", len_outputs)

	# Remove nan values
	for x in (inputs_masked):
		x[np.isnan(x)] = 0
		x[np.isinf(x)] = 0

	print("Input scaling: ", input_scale)
	print("Output scale min: ", output_scale_min)
	print("Output scale max: ", output_scale_max)

	# Rescale inputs masked
	eps = 1e-6
	for i, k in enumerate(input_scale):
		if k > eps:
			inputs_masked[:,i] /= k
		else:
			print("Non rescale because or zero values")

	# Load model       
	model = keras.models.load_model(model_path)
	batch_size = 1024

	# Apply model prediction to inputs masked
	print("prediction ...")
	y_predict = model.predict(inputs_masked, batch_size=batch_size)
	print("done")

	# Create output parameter folder
	try:
		if not os.path.exists(output_parameter_folder):
			os.makedirs(output_parameter_folder)
	except OSError:
		print ('Error: Creating directory. ' + output_parameter_folder)

	# Reshape output and save corresponding nifti 
	
	for k in range(0, len_outputs):
		output_predict_list = (y_predict[:,k]  * (output_scale_max[k] - output_scale_min[k])) + output_scale_min[k]
		
		output_predict = np.zeros(image_length)
		print("output_predict shape: ", output_predict.shape)

		output_predict[mask] = output_predict_list
		
		output_predict = output_predict.reshape(mask_shape[0], mask_shape[1], mask_shape[2])
		print("output_predict shape: ", output_predict.shape)
	 
		output_nii_name = base_signal_name + '_' + out_labels[k] + '_' + output_suffix + '.nii.gz'
		output_nii_path = os.path.join(output_parameter_folder, output_nii_name)	
		
		new_nii = nib.Nifti1Image(output_predict, affine, header)
		new_nii.to_filename(output_nii_path)

