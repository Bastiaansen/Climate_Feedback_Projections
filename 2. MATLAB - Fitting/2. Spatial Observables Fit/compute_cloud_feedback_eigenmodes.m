function compute_cloud_feedback_eigenmodes(Modes, save_file_name, tau, folder_path)
%compute_cloud_feedback_eigenmodes: perform linear regression to obtain the
%eigenmode decomposition for cloud feedback.
%   A loop over all lat and lon combinations is used and time series for
%   each of these location is multiplied against the inverse of the modes
%   (used as input) to obtain linear regression per location.

%% Obtain the relevant field files

field_name = '__xarray_dataarray_variable__';

% Planck
[dLW_PL, lat, lon, year] = loadData([folder_path 'dLW_PL_field.nc'], field_name);
[dLW_CS_PL, ~, ~, ~] = loadData([folder_path 'dLW_PL_field.nc'], field_name);

% Lapse Rate
[dLW_LR, ~, ~, ~] = loadData([folder_path 'dLW_lr_field.nc'], field_name);
[dLW_CS_LR, ~, ~, ~] = loadData([folder_path 'dLW_lr_field.nc'], field_name);

% Surface Albedo
[dSW_SA, ~, ~, ~] = loadData([folder_path 'dSW_alb_field.nc'], field_name);
[dSW_CS_SA, ~, ~, ~] = loadData([folder_path 'dSW_alb_field.nc'], field_name);

% Water Vapour
[dLW_WV, ~, ~, ~] = loadData([folder_path 'dLW_q_field.nc'], field_name);
[dLW_CS_WV, ~, ~, ~] = loadData([folder_path 'dLW_CS_q_field.nc'], field_name);

[dSW_WV, ~, ~, ~] = loadData([folder_path 'dSW_q_field.nc'], field_name);
[dSW_CS_WV, ~, ~, ~] = loadData([folder_path 'dSW_CS_q_field.nc'], field_name);

% Imbalance
[dN, ~, ~, ~] = loadData([folder_path 'dIMB_field.nc'], field_name);
[dN_CS, ~, ~, ~] = loadData([folder_path 'dIMB_CS_field.nc'], field_name);

%% Construct the cloud contribution
% CLOUD FB + (F - F_cs) = (dN - dN_cs) + sum( dFB_cs - dFB )

dCLOUD_PLUS_FORCING = (dLW_CS_PL - dLW_PL) + (dLW_CS_LR - dLW_LR) + ...
    (dSW_CS_SA - dSW_SA) + (dLW_CS_WV - dLW_WV) + (dSW_CS_WV - dSW_WV) + ...
    (dN - dN_CS);

%% Prepare for the regression
% Now also need a constant because of the (constant) forcing difference
% included in the time series

var_reshape = reshape(dCLOUD_PLUS_FORCING, length(lon)*length(lat), [], 1);

X = [ones(1,length(year)); Modes];

%% Regression
% In matlab this is just a single command

coeff_reshape = var_reshape / X;


%% Reshape coefficients
coeff = reshape(coeff_reshape, length(lon), length(lat), []);

cloud_modes = coeff(:,:,2:end); % First component corresponds to constant forcing difference
cloud_masking_of_forcing = coeff(:,:,1);

%% Save the coefficients/the eigenmode field
save_eigenmodes_var(save_file_name, cloud_modes, lon, lat, tau);

%% Save the cloud masking of forcing

save_file_name = char(save_file_name);
save_file_name = save_file_name(1:end-6);
save_cloud_masking([save_file_name(1:end-3) + "_masking_forcing.nc"], cloud_masking_of_forcing, lon, lat, 'Cloud_Masking')


