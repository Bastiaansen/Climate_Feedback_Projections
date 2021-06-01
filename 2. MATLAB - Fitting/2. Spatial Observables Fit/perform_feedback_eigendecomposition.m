%% Perform_feedback_eigendecomposition
% This script computes the different eigenmodes for the various climate
% feedback (and for surface temperature and radiative imbalance). For this,
% hard-coded results of eigenmode timescales tau needs to be put in below.
% Furthermore, in the current folder, 3D data needs to be put in that has
% latitude, longitude and temporal information of the feedback that needs
% to be split into its eigenmodes. The method here uses a linear regression
% for each (lat,lon) combination to fit the temporal evolution at that
% location to Y = c_1 e^(-t/tau_1) + ... + c_N e^(-t/tau_N). The
% coefficients c_1, ..., c_N are then combined back into a 2D field.

%% Clean-up matlabs workspace
close all
clear all

%% Path to the folder with the spatial feedback contributions
folder_path = '..\..\1. Python - CMIP Feedback Computations\CESM2 abrupt4xCO2\Data\fields\';

%% Prepare the planck feedback (by combining atmosphere & surface fields)
combine_Planck_surface_and_atmosphere_fields

%% Eigenmode timescales (hard-coded; result of different script)
tau = [ 4.5 ; 127 ; 889 ];

%% List of feedbacks/fields that need to be split into contributions
field_file_list = [ "dtas_field.nc"; 'dLW_PL_field.nc'; ...
    'dLW_lr_field.nc'; 'dSW_alb_field.nc'; 'dSW_q_field.nc'; 'dLW_q_field.nc'];
field_name_list = ["tas"; "__xarray_dataarray_variable__"; ...
    "__xarray_dataarray_variable__"; "__xarray_dataarray_variable__"; ...
    "__xarray_dataarray_variable__"; "__xarray_dataarray_variable__"];
eigenmode_file_name_suffic_list = [ "dtas"; 'Planck'; 'LR'; ...
    'alb'; 'q_SW'; 'q_LW'];

%% Construct X = [ e^(-t/tau) ] and compute its pseudo-inverse
% Then linear regression is only a multiplication of this pseudo-inverse to
% the imported data

% Import temporal data (years)
t = ncread([folder_path + field_file_list(1)], 'year');
t = double(t); % cast to double (since its int64 otherwise...)

% Eigenmodes have the form 1 - exp(-t/tau), ...
% except for the radiative imbalance where it has the form exp(-t/tau)

E = exp(- repmat(t', length(tau), 1) ./ tau);
X = 1 - E;

% Compute the pseudo-inverses
Einv = pinv(E);
Xinv = pinv(X);

%% Loop over all feedbacks and compute the eigenmode decomposition
for i = 1:length(field_file_list)
   compute_feedback_eigenmodes([folder_path + field_file_list(i)], field_name_list(i), Xinv, "Eigenmodes/" + string(eigenmode_file_name_suffic_list(i)) + "_modes.nc", tau);   
end

%% Also do that for the radiative imbalance (uses Einv instead of Xinv)
compute_IMB_feedback_eigenmodes([folder_path 'dIMB_field.nc'], "__xarray_dataarray_variable__", X, "Eigenmodes/dIMB_modes.nc", tau)

%% Finally, do the computations for the cloud feedback
% This needs its own function as cloud feedback is derived from the other
% feedbacks
compute_cloud_feedback_eigenmodes(X, "Eigenmodes/cloud_modes.nc", tau, folder_path)