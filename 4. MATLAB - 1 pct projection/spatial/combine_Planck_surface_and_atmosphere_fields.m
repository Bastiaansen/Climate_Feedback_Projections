%% Planck feedback (full-sky)

[dLW_PL_atm, lat, lon, years] = loadData([data_folder_path 'dLW_planck_atm_field.nc'], "__xarray_dataarray_variable__");
[dLW_PL_sur, ~, ~, ~] = loadData([data_folder_path 'dLW_ts_field.nc'], "__xarray_dataarray_variable__");

dLW_PL = dLW_PL_atm + dLW_PL_sur;

save_field([data_folder_path 'dLW_PL_field.nc'], dLW_PL, lon, lat, years)

%% Planck feedback (clear-sky)

[dLW_CS_PL_atm, ~, ~, ~] = loadData([data_folder_path 'dLW_CS_planck_atm_field.nc'], "__xarray_dataarray_variable__");
[dLW_CS_PL_sur, ~, ~, ~] = loadData([data_folder_path 'dLW_CS_ts_field.nc'], "__xarray_dataarray_variable__");

dLW_CS_PL = dLW_CS_PL_atm + dLW_CS_PL_sur;

save_field([data_folder_path 'dLW_CS_PL_field.nc'], dLW_CS_PL, lon, lat, years)


function save_field(file_name, var, lon, lat, year)

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

    %% Year
    nccreate(file_name, 'year', 'dimensions', {'year', length(year)});
    ncwriteatt(file_name, 'year', 'name', 'year');
    ncwrite(file_name, 'year', year);

    %% Field information
    nccreate(file_name, '__xarray_dataarray_variable__', 'dimensions', {'lon', length(lon), 'lat', length(lat), 'year', length(year)});
    ncwrite(file_name, '__xarray_dataarray_variable__' , var)
    
end