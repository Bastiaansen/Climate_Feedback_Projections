function [var, lat, lon, years] = loadData(file_name, field_name)

%% Load the data on lon, lat and the variable

lon = ncread(file_name, 'lon');
lat = ncread(file_name, 'lat');
years = ncread(file_name, 'year');
var = ncread(file_name, char(field_name));

%% Make sure var has dimensions lon, lat, year
var_meta = ncinfo(file_name, char(field_name));
var_dims = var_meta.Dimensions;

perms = zeros(1,length(var_dims));

for i = 1:length(var_dims)
    if(var_dims(i).Name == "lon")
        perms(1) = i;
    elseif(var_dims(i).Name == "lat")
        perms(2) = i;
    elseif(var_dims(i).Name == "year")
        perms(3) = i;
    elseif(var_dims(i).Name == "member_id")
        perms(4) = i;
    end
end
var = permute(var, perms);
var = squeeze(var);
    


end

