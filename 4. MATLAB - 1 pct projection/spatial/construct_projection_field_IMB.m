function construct_projection_field_IMB(data_file, field_name, mode_file, forcing_file_name, save_file_name)


%% Factor to transform 4xCO2 to 1pctCO2 experiment
beta_factor = log(1.01)/log(4);

%% Load the computed spatial modes
[modes, lat, lon, tau] = loadModes(mode_file);
modes_reshape = reshape(modes, length(lon)*length(lat),[], 1);

F = loadInitialForcing(forcing_file_name); % Forcing per
F_reshape = reshape(F, length(lon)*length(lat), [], 1);

%% Load the computed time series from the model output
[var, lat, lon, t] = loadData(data_file, field_name);
var_reshape = reshape(var, length(lon)*length(lat),[],1);
t = double(t)';


%% Compute projection time series
temporal_modes = repmat(t,length(tau),1) + tau .* ( exp(- repmat(t,length(tau),1) ./ abs(tau)) - 1);

proj_var_reshape = F_reshape * beta_factor * t + beta_factor * modes_reshape * temporal_modes;

error_reshape = proj_var_reshape - var_reshape;

%% Reshape back

proj_var = reshape(proj_var_reshape, length(lon), length(lat), []);
error = reshape(error_reshape, length(lon), length(lat), []);

%% Save the projected time series and the errors
save_projection(save_file_name, proj_var, error, lon, lat, t);


end

