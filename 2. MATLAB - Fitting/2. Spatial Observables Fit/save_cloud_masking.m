function save_cloud_masking(file_name, cloud_masking, lon, lat, var_name)
%save_cloud_maksing: saves the found cloud masking of forcing

file_name = char(file_name)
%% Delete existing file if it exists
if isfile(file_name)
   delete(file_name) 
end

%% Longitude
nccreate(file_name, 'lon', 'dimensions', {'lon', length(lon)});
ncwriteatt(file_name, 'lon', 'standard_name', 'longitude');
ncwriteatt(file_name, 'lon', 'units', 'degrees_east');
ncwriteatt(file_name, 'lon', 'valid_max', 360);
ncwriteatt(file_name, 'lon', 'valid_min', 0);
ncwrite(file_name, 'lon', lon);

%% Latitude
nccreate(file_name, 'lat', 'dimensions', {'lat', length(lat)});
ncwriteatt(file_name, 'lat', 'standard_name', 'latitude');
ncwriteatt(file_name, 'lat', 'units', 'degrees_north');
ncwriteatt(file_name, 'lat', 'valid_max', 90);
ncwriteatt(file_name, 'lat', 'valid_min', -90);
ncwrite(file_name, 'lat', lat);


%% Eigenmodes
nccreate(file_name, var_name, 'dimensions', {'lon', length(lon), 'lat', length(lat)});
ncwrite(file_name, var_name, cloud_masking)


end

