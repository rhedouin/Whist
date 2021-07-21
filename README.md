# Whist
White matter simulation toolbox

This toolbox is related to a paper available on Neuroimage
https://www.sciencedirect.com/science/article/pii/S1053811921004158

This matlab toolbox has 2 main purposes 

(1) Create 2D WM models with chosen microstructure parameters using realistic axon shapes. 

The corresponding function are placed in the createWMModel folder. \
See createOneWMModelExample.m for an example of usage \
See complete_WM_Model_creation.mp4 for a video illustration

(2) Simulate the field perturbation and the corresponding multi GRE signal from a 2D or 3D WM model. 

The corresponding function are placed in the signalSimulation folder. \
See simulateFieldPerturbationExample.m for an example of usage \
See simulate_field_perturbation_and_GRE_signal.mp4 for an illustration

## Terms of use 
If you have any question or you would like to provide suggestion to improve this toolbox/report bug(s) please feel free to contact me renaud.hedouin@gmail.com \
All the codes and methods developed in Whist are under MIT license. \
The axon packing part is based on a toolbox developed by Tom Mingasson for hollow cylinder models, thanks to him \
https://github.com/neuropoly/axonpacking

Renaud
