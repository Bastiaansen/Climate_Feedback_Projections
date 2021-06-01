% This script can be used to create spatial projections for the 1pctCO2
% scenario. It takes care of
% (1) loading in data derived from the actual GCM
% (2) loading in the fits from the abrupt4xCO2 experiment
% (3) creating projections for the various spatial observables in the
% 1pctCO2 experiment
% (4) saving the computed projections

%% clean slate

clear all
close all

%% Path to the folder with the spatial feedback contributions
data_folder_path = '..\..\1. Python - CMIP Feedback Computations\CESM2 1pct CO2\Data\fields\';

%% Prepare the planck feedback (by combining atmosphere & surface fields)
combine_Planck_surface_and_atmosphere_fields

%% List of feedbacks/fields
field_file_list = [ "dtas_field.nc"; 'dLW_PL_field.nc'; ...
    'dLW_lr_field.nc'; 'dSW_alb_field.nc'; 'dSW_q_field.nc'; 'dLW_q_field.nc'];
field_name_list = ["tas"; "__xarray_dataarray_variable__"; ...
    "__xarray_dataarray_variable__"; "__xarray_dataarray_variable__"; ...
    "__xarray_dataarray_variable__"; "__xarray_dataarray_variable__"];
eigenmode_file_name_suffic_list = [ "dtas"; 'Planck'; 'LR'; ...
    'alb'; 'q_SW'; 'q_LW'];

%% Path to the folder with the eigenmodes (as computed from abrupt4xCO2)

mode_folder_path = '..\..\2. MATLAB - Fitting\2. Spatial Observables Fit\Eigenmodes\';

%% Path to the folder in which the projections need to be stored

save_folder_path = 'spatial_projection\';

%% Construct the projected states for each spatial observable

for i = 1:length(field_file_list)
    construct_projection_field([data_folder_path + field_file_list(i)], ...
        field_name_list(i), ...
        [mode_folder_path + eigenmode_file_name_suffic_list(i) + "_modes.nc"], ...
        [save_folder_path + eigenmode_file_name_suffic_list(i) + "_projection.nc"])
end

%% Also do that for the radiative imbalance (different because of forcing)
construct_projection_field_IMB([data_folder_path + "dIMB_field.nc"], ...
    "__xarray_dataarray_variable__", [mode_folder_path + "dIMB_modes.nc"], ...
    [mode_folder_path + "dIMB_initial_forcing.nc"], ...
    [save_folder_path + "dIMB_projection.nc"]);

%% Finally, also for the cloud feedback
construct_projection_field_CLOUD(data_folder_path, mode_folder_path, ...
    save_folder_path);