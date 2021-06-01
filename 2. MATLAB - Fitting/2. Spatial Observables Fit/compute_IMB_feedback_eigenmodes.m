function compute_IMB_feedback_eigenmodes(file_name, field_name, Modes, save_file_name, tau)
%compute_feedback_eigenmodes: perform linear regression to obtain the
%eigenmode decomposition of a (feedback) field.
%   A loop over all lat and lon combinations is used and time series for
%   each of these location is multiplied against the inverse of the modes
%   (used as input) to obtain linear regression per location.

%% Load the data on lon, lat and the variable

[var, lat, lon, ~] = loadData(file_name, field_name);

%% Reshape variable
% The variable is a 3D array with (lon, lat, time) dimensions, but we need
% to reshape it to a 2D array with (location, time) dimensions
var_reshape = reshape(var, length(lon)*length(lat), [], 1);

%% Prepare for the regression
% Now also need a constant because of the (constant) forcing difference
% included in the time series

X = [ones(1,length(Modes(1,:))); Modes];

%% Perform linear regression
% Since we have var_reshape = [ y(t1), ..., y(tM)] and Xinv this is now a
% matrix multiplication

coeff_reshape = var_reshape / X;

%% Reshape coefficients
% We now a 2D array for the regression coefficients with dimensions
% (location, coefficients), but we want to go back to a 3D array that has
% dimensions (lon, lat, coefficients), so we need to do another reshape
coeff = reshape(coeff_reshape, length(lon), length(lat), []);

IMB_modes = coeff(:,:,2:end); % First component is (initial) forcing
IMB_forcing = coeff(:,:,1);

%% Save the coefficients/the eigenmode field
save_eigenmodes_var(save_file_name, IMB_modes, lon, lat, tau)

%% Save the initial forcing
save_file_name = char(save_file_name);
save_file_name = save_file_name(1:end-6);
save_cloud_masking([save_file_name(1:end-3) + "_initial_forcing.nc"], IMB_forcing, lon, lat, 'Initial_Forcing')



end

