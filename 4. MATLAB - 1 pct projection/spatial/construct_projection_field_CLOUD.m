function construct_projection_field_CLOUD(data_folder, mode_folder, save_folder)

%% Obtain the relevant field files

field_name = '__xarray_dataarray_variable__';

% Planck
[dLW_PL, lat, lon, t] = loadData([data_folder 'dLW_PL_field.nc'], field_name);
[dLW_CS_PL, ~, ~, ~] = loadData([data_folder 'dLW_PL_field.nc'], field_name);

% Lapse Rate
[dLW_LR, ~, ~, ~] = loadData([data_folder 'dLW_lr_field.nc'], field_name);
[dLW_CS_LR, ~, ~, ~] = loadData([data_folder 'dLW_lr_field.nc'], field_name);

% Surface Albedo
[dSW_SA, ~, ~, ~] = loadData([data_folder 'dSW_alb_field.nc'], field_name);
[dSW_CS_SA, ~, ~, ~] = loadData([data_folder 'dSW_alb_field.nc'], field_name);

% Water Vapour
[dLW_WV, ~, ~, ~] = loadData([data_folder 'dLW_q_field.nc'], field_name);
[dLW_CS_WV, ~, ~, ~] = loadData([data_folder 'dLW_CS_q_field.nc'], field_name);

[dSW_WV, ~, ~, ~] = loadData([data_folder 'dSW_q_field.nc'], field_name);
[dSW_CS_WV, ~, ~, ~] = loadData([data_folder 'dSW_CS_q_field.nc'], field_name);

% Imbalance
[dN, ~, ~, ~] = loadData([data_folder 'dIMB_field.nc'], field_name);
[dN_CS, ~, ~, ~] = loadData([data_folder 'dIMB_CS_field.nc'], field_name);


%% Construct the cloud contribution

dCLOUD_PLUS_FORCING = (dLW_CS_PL - dLW_PL) + (dLW_CS_LR - dLW_LR) + ...
    (dSW_CS_SA - dSW_SA) + (dLW_CS_WV - dLW_WV) + (dSW_CS_WV - dSW_WV) + ...
    (dN - dN_CS);

var_reshape = reshape(dCLOUD_PLUS_FORCING, length(lon)*length(lat), [], 1);
t = double(t)';



%% Factor to transform 4xCO2 to 1pctCO2 experiment
beta_factor = log(1.01)/log(4);

%% Load the computed spatial modes
[modes, lat, lon, tau] = loadModes([mode_folder + "cloud_modes.nc"]);
modes_reshape = reshape(modes, length(lon)*length(lat),[], 1);

[F_CLOUD, lat, lon] = loadCloudMasking([mode_folder + "cloud_masking_forcing.nc"]);
F_CLOUD_reshape = reshape(F_CLOUD, length(lon)*length(lat), [], 1);


%% Compute projection time series
temporal_modes = repmat(t,length(tau),1) + tau .* ( exp(- repmat(t,length(tau),1) ./ abs(tau)) - 1);

proj_var_reshape = beta_factor * modes_reshape * temporal_modes;

var_minus_forcing_reshape = var_reshape - F_CLOUD_reshape .* beta_factor .* ones(length(F_CLOUD_reshape),length(t')).*t;

error_reshape = proj_var_reshape - var_minus_forcing_reshape;


%% Reshape back

proj_var = reshape(proj_var_reshape, length(lon), length(lat), []);
error = reshape(error_reshape, length(lon), length(lat), []);
var_minus_forcing = reshape(var_minus_forcing_reshape, length(lon), length(lat), []);

%% Save the projected time series and the errors AND the data-cloud_masked_forcing
save_projection_CLOUD([save_folder + "cloud_projection.nc"], proj_var, error, var_minus_forcing, lon, lat, t);

end

