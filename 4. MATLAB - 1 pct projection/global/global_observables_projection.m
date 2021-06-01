% This script can be used to create projections for the global observables
% in the 1pctCO2 experiment. It takes care of the following:
% (1) load the actual data that is computed from the GCM run
% (2) load the fit made from the abrupt4xCO2 experiment
% (3) create projection for the 1pctCO2 experiment
% (4) make some plots, showing evolution over time

%% clean slate

clear all
close all

%% Load-in data
var_names = '__xarray_dataarray_variable__';
folder_path = '..\..\1. Python - CMIP Feedback Computations\CESM2 1pct CO2\Data\global\';

GMST = ncread([folder_path  'dGMST.nc'], var_names);
IMB = ncread([folder_path 'dIMB.nc'], var_names);

dLW_LR = ncread([folder_path 'dLW_lr.nc'], var_names);
dLW_PL = ncread([folder_path 'dLW_planck.nc'], var_names);
dSW_AL = ncread([folder_path 'dSW_alb.nc'], var_names);
dLW_WV = ncread([folder_path 'dLW_q.nc'], var_names);
dSW_WV = ncread([folder_path 'dSW_q.nc'], var_names);

% Clear sky contributions to radiative imbalance of feedbacks
dLW_CS_PL = ncread([folder_path 'dLW_CS_planck.nc'], var_names);
dLW_CS_LR = ncread([folder_path 'dLW_CS_lr.nc'], var_names);
dSW_CS_AL = ncread([folder_path 'dSW_CS_alb.nc'], var_names);
dLW_CS_WV = ncread([folder_path 'dLW_CS_q.nc'], var_names);
dSW_CS_WV = ncread([folder_path 'dSW_CS_q.nc'], var_names);

IMB_CS = ncread([folder_path 'dIMB_CS.nc'], var_names);

dLW_ts = ncread([folder_path 'dLW_ts.nc'], var_names);
dLW_ta = ncread([folder_path 'dLW_ta.nc'], var_names);

t = double(ncread([folder_path 'dGMST.nc'], 'year'));

%% Compute cloud feedback plus forcing difference
dCLOUD_PLUS_FORCING = (dLW_CS_PL - dLW_PL) + (dLW_CS_LR - dLW_LR) + ...
    (dSW_CS_AL - dSW_AL) + (dLW_CS_WV - dLW_WV) + (dSW_CS_WV - dSW_WV) + ...
    (IMB - IMB_CS);

%% Load-in fitted parameters
fit_data_folder = '../../2. MATLAB - Fitting/1. Global Observables Fit/fit_data';
load([fit_data_folder '/best_fits.mat'])
% We use the M = 3 scenario here
M = 3;
x = x_fits{3};
conf = confidence_intervals{3};

%% Obtain (projected) time series for the 1pct scenario from fitted params
t_fit = [0:1:200]';
Y = exponent_fit_function_1pct(x,t_fit');

%% PLOTTING of projections plus data

plot_data_plus_projection

%% PLOTTING of error between data and projection

plot_projection_errors