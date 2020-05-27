#!/usr/bin/env python
# coding: utf-8


def quick_sequential(
		input_shape,
		output_dim,
		layers=(
				(lambda i, o: int(2 * i * o), 'tanh', 0.4, False, False),
				(lambda i, o: int(1.5 * i * o), 'tanh', 0.2, True, False),
				(lambda i, o: int(1.25 * i * o), 'tanh', 0.1, False, False),
				(None, 'linear', 0.0, False, False))):

	try:
		iter(input_shape)
	except TypeError:
		input_shape = (input_shape,)
	input_dim = np.sum(np.array(input_shape))

	model = keras.models.Sequential()

	# input layer
	model.add(keras.layers.Dense(input_dim, input_shape=input_shape))
	model.add(keras.layers.normalization.BatchNormalization(
			 input_shape=input_shape))

	name = 'mlp__{input_dim}'.format_map(locals())

	# hidden and output layer
	for layer in layers:
		dim_func, activation, dropout, normalization, bias = layer
		flags = ((
			('N', normalization),
			('B', bias)))
		out_dim = dim_func(input_dim, output_dim) if dim_func else output_dim
		model.add(keras.layers.Dense(
			# out_dim, kernel_initializer='lecun_uniform', use_bias=bias, kernel_regularizer=l2(0.001)))
		out_dim, kernel_initializer='lecun_uniform', use_bias=bias))

		if activation:
			model.add(keras.layers.Activation(activation))
		if dropout:
			model.add(keras.layers.Dropout(dropout))
		if normalization:
			model.add(keras.layers.normalization.BatchNormalization())

		name += '_{out_dim}-{activation}-{dropout}'.format_map(locals())
		if any([v for k, v in flags]):
			name += '_'
			for k, v in flags:
				if v:
					name += k

	sgd = keras.optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
	loss='mean_absolute_error'
						  
	name += '__{sgd}__{loss}'.format_map(locals())
	name += '__{batch_size}__{epochs}.{keras_backend}.hdf5'
	
	model.compile(optimizer=sgd, loss=loss)   
	
	return model, name

# Load libraries
import os

import numpy as np
import scipy as sp

import scipy.io

import keras
import tensorflow

from keras.regularizers import l1, l2

import h5py
import matplotlib.pyplot as plt

from shutil import copyfile
import time

import json

# Set parameters
load_model = 0
if load_model == 1:
	input_nb_epoch =15 # if load model
else:
	input_nb_epoch = 0
	
nb_epochs_to_perform = 40 # to perform on the load model or the generated model

batch_size = 1024

######## Need to be fill
# Dictionary details
nb_TE = 18
length_one_orientation = 2*nb_TE + 1
nb_orientations = 10
nb_replic = 8
noise = '1'
# Name project / experiment
base_folder = '/project/3015069.04/'
experiment_name = 'BS-3'
########

print("noise level : " + noise)
	
# Name Input/Output data

dict_folder  = os.path.join(base_folder, 'dictionaries/multi_orientations/', experiment_name)
input_base_signal_name = 'SignalWithNoise' + noise + '_' + str(nb_replic) + 'rep_' +  str(nb_orientations) + 'rotations_' +  str(nb_TE) + 'TE_' + experiment_name + '_polyfit_cartesian_with_theta'
signal_filepath = os.path.join(dict_folder, input_base_signal_name + '.h5py')

output_base_signal_name = input_base_signal_name 

base_dp_folder = os.path.join(base_folder, 'deep_learning/multi_orientations/', experiment_name)
dp_folder = os.path.join(base_dp_folder, output_base_signal_name)

try:
   if not os.path.exists(dp_folder):
	   os.makedirs(dp_folder)
except OSError:
   print ('Error: Creating directory. ' + dp_folder)
	   
input_scale_name = output_base_signal_name + '_input_scale.npy'
input_scale_path = os.path.join(dp_folder,input_scale_name)

output_scale_min_path = os.path.join(dp_folder,output_base_signal_name + '_output_scale_min.npy')
output_scale_max_path = os.path.join(dp_folder,output_base_signal_name + '_output_scale_max.npy')

# Load Input/Output data
print('Load data ...')

mat = h5py.File(signal_filepath)

labels = 'SignalValues', 'FVFValues', 'gRatioValues', 'xiMyelinValues', 'T2IntraExtraAxonalValues', 'T2MyelinValues', 'weightValues'
out_labels = 'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'

composite_at_begin_labels = 'SignalValues',
dico_shape = mat[labels[1]].shape

data = np.concatenate(
	[np.array(mat[label]) for label in composite_at_begin_labels] + 
	[np.array(mat[label]).reshape(dico_shape + (1,)) for label in labels
	 if label not in composite_at_begin_labels],-1) 

data = np.float32(data)

len_inputs = nb_orientations * length_one_orientation 
len_outputs = data.shape[-1] - len_inputs
print("Len Inputs: ", len_inputs)
print("Len Outputs: ", len_outputs)

# Data preprocessing
print('Preprocess data ...')
print("data shape : ", data.shape)
data = data.reshape(-1,data.shape[-1])
print("data shape : ", data.shape)

inputs = data[:,:len_inputs].copy()
outputs = data[:,len_inputs:].copy()

print("inputs shape : ", inputs.shape)

del data

mask = np.ones(dico_shape)
mask = mask.reshape(-1)
mask = mask.astype(bool)

print("dico_shape: ", dico_shape)
print("mask_shape: ", mask.shape)

for x in (inputs, outputs):
	x[np.isnan(x)] = 0
	x[np.isinf(x)] = 0
# scale data
tail = 0.05
eps = 1e-6

input_scaling = []
for i in range(inputs.shape[1]):
	k = np.max(
		(np.abs(np.percentile(inputs[:,i], tail * 100)),
		 np.abs(np.percentile(inputs[:,i], (1 - tail) * 100))))
	input_scaling.append(k)
	if k > eps:
		inputs[:,i] /= k
	else:
		print("avoid division by 0")
		
	
print("Save input scale: ", input_scaling)
np.save(input_scale_path, input_scaling)

output_scaling_max = []
output_scaling_min = []
output_scaling_diff = []

for i in range(outputs.shape[1]):
	k_max = np.max(outputs[:,i])
	k_min = np.min(outputs[:,i])
	k_diff = k_max - k_min

	output_scaling_max.append(k_max)
	output_scaling_min.append(k_min)
	output_scaling_diff.append(k_diff)
	
	outputs[:,i] = (outputs[:,i] - k_min) / (k_diff)
	
print("Save output scale diff: ", output_scaling_diff)
np.save(output_scale_max_path, output_scaling_max)
np.save(output_scale_min_path, output_scaling_min)

# Mask data
x_train = inputs[mask, :]
y_train = outputs[mask, :]

print("x_train shape : ", x_train.shape)
print("y_train shape : ", y_train.shape)

print('done')
	
# Create or load model
if load_model == 0:
	print('create model ...')	   
	model, ann_name = quick_sequential(len_inputs, len_outputs)
	keras_backend = keras.backend.backend()
	ann_filepath = os.path.join(dp_folder, 'epoch_' + str( nb_epochs_to_perform) + '.h5py')
   
else:
	model_name = 'epoch_' + str(input_nb_epoch) + '.h5py'
	input_model_path = os.path.join(dp_folder, model_name)
	print('load model ', input_model_path)
	model = keras.models.load_model(input_model_path)
	ann_filepath = os.path.join(dp_folder, 'epoch_' + str(input_nb_epoch + nb_epochs_to_perform) + '.h5py')
  
print('done')


print('ANN training...')
	
# Deep learning training
print('ANN: {}'.format(ann_filepath))

callbacks= [
	keras.callbacks.EarlyStopping(monitor='loss',patience=5),
	keras.callbacks.TensorBoard(log_dir=os.path.join(base_folder, 'logs')),
	keras.callbacks.ModelCheckpoint(ann_filepath + '.check', save_best_only=True, save_weights_only=False)]

history_callback = model.fit(
	x_train, y_train, batch_size=batch_size, epochs=nb_epochs_to_perform,
	validation_split=1/nb_replic, callbacks=callbacks, shuffle=True)

model.save(ann_filepath)

# Save the history under the form of a json file
history_path = os.path.join(dp_folder, 'total_history.json')
history_dict = history_callback.history
json.dump(history_dict, open(history_path, 'w'))
		
