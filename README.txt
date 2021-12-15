Whist: WHIte matter Stimulation Toolbox

This toolbox was developped to create WM model and use it to obtain microstructural information.
Refer to this paper:

Decoding the microstructural properties of white matter using realistic
models, NeuroImage, Hédouin et al.

This toolbox is associated with a data sharing collection to allow to reproduce most of the results of the article
https://data.donders.ru.nl/collections/di/dccn/DSC_3015069.04_445?0

This toolbox has the following functionalities:

create WM Models (1)
simulate the corresponding ME-GRE (2)
create a dictionary from these signals (3)
train a deep learning with the dictionary (4)

post-process data (5)
recover brain parameter maps (6)


(1) Create 2D WM models with chosen microstructure parameters using realistic axon shapes.

See createOneWMModelExample.m for an example of usage
Easy to parallelize if a large number of WM models are needed
See complete_WM_Model_creation.mp4 for a video illustration

Another option is now available to create one axon with myelin water layers (see WMModelMyelinWithLayer). Several functions directly create the figures presented in Appendix B of the paper.

(2) Simulate the field perturbation and the corresponding ME-GRE signal from a 2D or 3D WM model.

See simulateFieldPerturbationExample.m for an example of usage
See simulate_field_perturbation_and_GRE_signal.mp4 for an illustration

(3) Simulate the ME-GRE signal for many WM models and microstructure parameters range

everything is parallelized
createDictionaryCompletePipelineWithLorentzianCorrection.m directly creates the dictionary from the WM models following several steps:
computeCompartmentSignalFromModels creates the field perturbation and the corresponding ME-GRE signals for different WM models parameters (FVF, g-ratio) for each compartment (intra-axonal, extra-axonal, myelin) 
computeTotalSignalDictionaryPart creates the total signal adding the T2 weighting, the compartment weights and normalize the signal
concatenateDictionaryParts concatenates the previous results into an entire dictionary

(4) Deep_learning training
Keras_train_WMModel_BS3_all_orientations_with_regularization_new_scaling.py within deep_learning
To be used on a GPU, just need to fill the dictionary details

(5) exVivoProcessing
see transformAndConcatenateBS3DataMain.m (parallalized along flip angles)
Processed the data from acquisition, take the magn and unwrapped phase all registered in the same space (same rotations should be used than with the dictionary) and transform it into a ME-GRE signal normalized concatenate along all orientations and include the theta fiber orientation.

(6) Deep_learning application to processed data
recover_parameter_BS3_all_orientations_polyfit_cartesian_new_scaling.py
To be used on GPU, create brain parameter maps from the trained deep learning (4) and the concatenate signals (5)

If you have any question or you would like to provide suggestion to improve this toolbox/report bug(s) please feel free to contact me renaud.hedouin@gmail.com (Renaud Hédouin)


