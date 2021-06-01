function save_projection(file_name, proj_var, error, lon, lat, year)
% save_projection: saves the found spatial projection and the resulting
% discrepancy/error between the projection and the feedback contributions
% computed directly from the model output.

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

%% time
nccreate(file_name, 'year', 'dimensions', {'year', length(year)});
ncwriteatt(file_name, 'year', 'name', 'year')
ncwrite(file_name, 'year', year)

%% Projection
nccreate(file_name, 'Projection', 'dimensions', {'lon', length(lon), 'lat', length(lat), 'year', length(year)});
ncwrite(file_name, 'Projection', proj_var)

%% Error
nccreate(file_name, 'Error', 'dimensions', {'lon', length(lon), 'lat', length(lat), 'year', length(year)});
ncwrite(file_name, 'Error', error)

end
