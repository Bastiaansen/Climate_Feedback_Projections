function save_eigenmodes_var(file_name, eigenmodes, lon, lat, tau)
%save_eigenmodes_var: saves the found eigenmodes in the cells of the
%variable eigenmodes to a nc file with seperate entry for each of the
%eigenmodes. Eigenmodes are 3D arrays with first two dimensions
%corresponding to longitude and latitude; the third dimensions corresponds
%to the different modes.

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

%% Taus
nccreate(file_name, 'tau', 'dimensions', {'tau', length(tau)});
ncwriteatt(file_name, 'tau', 'name', 'Eigenmode timescale')
ncwrite(file_name, 'tau', tau)

%% Eigenmodes
nccreate(file_name, 'Eigenmodes', 'dimensions', {'lon', length(lon), 'lat', length(lat), 'tau', length(tau)});
ncwrite(file_name, 'Eigenmodes', eigenmodes)


end

