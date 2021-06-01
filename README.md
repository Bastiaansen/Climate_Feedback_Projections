# Projections of the Transient State-Dependency of Climate Feedbacks


## Instructions and Explanations of codes used

Explanation is per map.



### 1. Python

For computations of observables from GCM data (2D fields and globally averaged), which gets saved in nc files. These require radiative kernels. For this study, we have used CAM5 kernels by Pendergrass ..... The Jupyter notebooks should be run to create these data files. Most of the actual code is contained in seperate python files.

There are two submaps for the different forcing scenarios:
	* CESM2 abrupt4xCO2
	* CESM2 1pctCO2
	
	
### 2. MATLAB - Fitting

Used to fit previously obtained observable data on feedback contributions, net radiative imbalance and temperature from the abrupt4xCO2 experiment to the response function (a sum of exponential decaying functions) -- see the paper for more details.

Submaps:
1. **global observables fit**

	Used for fits of the globally averaged observables.

	To create the fits, use the script 'compute_feedbacks.m'. This script automatically finds the computed observable files if left in their original folder. The code computes cloud feedback contributions and tries values M=1 to M=5 for the amount of (Koopman) modes, create some statistics of these tests and saves some fits (in the folder 'fit-data'). It also creates plots of the best found fits for different values for M (i.e. the amount of modes). Finally, M = 3 is taken to create the plots in the paper and to compute feedback strengths per mode.

	This script relies on the other MATLAB files in the folder:
	* exponent_fit_function.m: contains the function to be fitted to the data
	* exponent_fit_function_der.m: contains the derivative of the function to be fitted (used for computations of instantaneous feedback strengths)
	* compute_feedbacks_for_fit.m: computes feedback strengths per mode and equilibrium feedback strengths from fitted parameters (contained in the variable x). Also propogates uncertainties.
	* propogate_confidenceintervals_sum: propogation of uncertainties for a sum x_1 + ... + x_n.
	* propogate_confidenceintervals_quotient.m: propogation of uncertainties for fractions p/q.
	* append_plot_with_fits.m: add fitted function to plots of observables (cloud feedback included masking of forcing).
	* append_plot_with_fit_clouds.m: add fitted function to plots of observables (but cloud feedback is seperated from masking of forcing).
		
2. **spatial observables fit**

	For fits of the spatial modes. The requires timescales tau_m obtained from nonlinear regression on globally averaged observables.

	To obtain the spatial fits, the script 'perform_feedback_eigendecomposition.m' should be run. It automatically loads the 2D fields if left in original folder, prepares those files and computes 2D modes from linear regression (coded here via multiplication with Moore-Penrose pseudoinverse matrices). Computed 2D modes are saved in the folder 'Eigenmodes'.

	This script depends on other MATLAB files in the folder:
	* loadData.m: load the 2D fields and prepare them for analysis.
	* combine_Planck_surface_and_atmosphere_fields.m: combines the fields that constitute the Planck feedback.
	* compute_feedback_eigenmodes.m: compute and save 2D modes for generic feedback.
	* compute_IMB_feedback_eigenmodes.m: compute and save 2D modes for radiative imbalance (and initial radiative forcing).
	* compute_cloud_feedback_eigenmodes.m: compute and save 2D modes for cloud feedback (and cloud masking of forcing).
	* save_eigenmodes_var.m: save the 2D modes found to a nc file.
	* save_cloud_masking.m: save the cloud masking found to a nc file.
		
3. **Analysis**

	For computation of (instantaneous) feedback strengths over time in the abrupt4xCO2 and the 1pctCO2 experiments, via the script 'Feedback_dissected.m'. It depends on some scripts in the folder '1. global observables fit' and should be run while in that folder.
		
### 3. Python Visualisations
For visualisations of the 2D modes, estimated equilibrium 2D spatial paterning in the abrupt4xCO2 epxeriment (both via the Jupyter notebook 'make_mode_plots.ipynb') and 2D projections compared to actual data for the 1pctCO2 experiment (via the Jupyter notebook 'compare 1pct_projections.ipynb').

### 4. MATLAB - 1pct projection
For creation of projections for the 1pctCO2 experiment, made from the fits of the abrupt4xCO2 experiment.

This contains two subfolders
1. **global**

	For projections of the globally averaged observables. The script 'global_observables_projection.m' takes care of (1) loading in the actual data for the 1pctCO2 experiment to compare projections with (files should be loaded automatically if left in place), (2) loading the fitted parameters, (3) creation of projections and (4) making some plots.

	The script relies on other MATLAB files in the folder:
	* exponent_fit_function_1pct.m: the response function for the 1pctCO2 experiment.
	* plot_data_plus_projection.m: plots data and projections.
	* add_LRT_projection.m: add linear response theory proections made with numeric (i.e. non-fitted) data Green functions (not shown in paper)
	* plot_projection_error.m: plot the mismatch of the projections (not shown in paper).

2. **spatial**

	For projections of the 2D fields of the observables. To run, use the script 'spatial_observables_projection.m' which loads in the data from the actual GCM and the fitted spatial modes. Also constructs spatial projections for the observables, which are saved to nc files in the folder 'spatial projection'.

	The script relies on other MATLAB files in the folder:
	* loadData.m: to load data derived from teh actual GCM 1pctCO2 experiment.
	* loadModes.m: to load the fitted spatial modes.
	* loadInitialForcing.m: to load the fitted effective initial radiative forcing.
	* loadCloudMasking.m: to load the fitted cloud masking of radiative forcing.
	* combine_Planck_surface_and_atmosphere_fields.m: combines the fields that constitute the Planck feedback.
	* construct_projection_field.m: construct 2D projection for generic climate feedbacks.
	* construct_projection_field_CLOUD.m: construct 2D projection for cloud feedback field, and compute (derived) actual cloud feedback.
	* construct_projection_field_IMB.m: construct projection for the net radiative imbalance.
	* save_projection.m: saves the projection to a nc file.
	* save_projection_cloud.m: saves the cloud feedback to a nc file.
		
### Python Packages

		
