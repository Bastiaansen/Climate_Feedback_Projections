function compute_feedback_eigenmodes(file_name, field_name, ModeInv, save_file_name, tau)
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

%% Perform linear regression
% Since we have var_reshape = [ y(t1), ..., y(tM)] and Xinv this is now a
% matrix multiplication

coeff_reshape = var_reshape * ModeInv;

%% Reshape coefficients
% We now a 2D array for the regression coefficients with dimensions
% (location, coefficients), but we want to go back to a 3D array that has
% dimensions (lon, lat, coefficients), so we need to do another reshape
coeff = reshape(coeff_reshape, length(lon), length(lat), []);

%% Save the coefficients/the eigenmode field
save_eigenmodes_var(save_file_name, coeff, lon, lat, tau)


end

