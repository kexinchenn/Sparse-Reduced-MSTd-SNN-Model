# Sparse, reduced representations in model of MSTd

This is the code accompanying the following publication:

> Chen, K., Beyeler, M., Krichmar, J. L. (2022) Cortical Motion Perception Emerges from Dimensionality Reduction with Evolved Spike-Timing Dependent Plasticity Rules. The Journal of Neuroscience (in Press). 


## Requiremments
- CARLsim5: https://github.com/UCI-CARL/CARLsim5
- ECJ: https://github.com/GMUEClab/ecj
- MATLAB 2017a

## Generate input stimuli
Stimuli in the training, validation, and test datasets can be generated with the MATLAB script `matlab_analysis_scripts/GenerateInputStim/generateInputStim.m`.

## Network Simulation
### evolve_MST_SNN_model_CARLsim
- Project code to evolve network parameters.
- Compile with `make`, and use `./launchCARLsimECJ.sh` to launch evolutionary runs. Network fitness values and evolved hyper-parameters are saved to out.stat.

### test_MST_SNN_model_CARLsim
- Project code that includes additional test trials.
- Pass in evolved parameters to simulate individual networks. Network activity and trial indices are saved in the `results` folder.

## Analysis scripts
Folder `matlab_analysis_scripts` contains MATLAB scripts for data analysis and visualization. These scripts were partially based on code developed for Beyeler et al. (2016).

## References
Beyeler, M., Dutt, N., and Krichmar, J.L. (2016). 3D Visual Response Properties of MSTd Emerge from an Efficient, Sparse Population Code. The Journal of Neuroscience 36, 8399-8415.