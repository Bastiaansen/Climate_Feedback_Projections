function [var, lat, lon] = loadInitialForcing(file_name)

%% Load the data on lon, lat and the variable

lon = ncread(file_name, 'lon');
lat = ncread(file_name, 'lat');
var = ncread(file_name, 'Initial_Forcing');

%% Make sure var has dimensions lon, lat, year
var_meta = ncinfo(file_name, 'Initial_Forcing');
var_dims = var_meta.Dimensions;

perms = zeros(1,length(var_dims));

for i = 1:length(var_dims)
    if(var_dims(i).Name == "lon")
        perms(1) = i;
    elseif(var_dims(i).Name == "lat")
        perms(2) = i;
    end
end
var = permute(var, perms);
var = squeeze(var);
    


end

