%% Adds projections made with (numerical) linear response theory
% That is, no fit is being made to the abrupt4xCO2 run.
% This shows the discrepancy between the computations with and without a
% fit.

%% Load-in data of abrupt4xCO2 experiment

var_names = '__xarray_dataarray_variable__';
folder_path = '..\..\1. Python - CMIP Feedback Computations\CESM2 abrupt4xCO2\Data\global\';

abr_GMST = ncread([folder_path  'dGMST.nc'], var_names);
abr_IMB = ncread([folder_path 'dIMB.nc'], var_names);

abr_dLW_LR = ncread([folder_path 'dLW_lr.nc'], var_names);
abr_dLW_PL = ncread([folder_path 'dLW_planck.nc'], var_names);
abr_dSW_AL = ncread([folder_path 'dSW_alb.nc'], var_names);
abr_dLW_WV = ncread([folder_path 'dLW_q.nc'], var_names);
abr_dSW_WV = ncread([folder_path 'dSW_q.nc'], var_names);

% Clear sky contributions to radiative imbalance of feedbacks
abr_dLW_CS_PL = ncread([folder_path 'dLW_CS_planck.nc'], var_names);
abr_dLW_CS_LR = ncread([folder_path 'dLW_CS_lr.nc'], var_names);
abr_dSW_CS_AL = ncread([folder_path 'dSW_CS_alb.nc'], var_names);
abr_dLW_CS_WV = ncread([folder_path 'dLW_CS_q.nc'], var_names);
abr_dSW_CS_WV = ncread([folder_path 'dSW_CS_q.nc'], var_names);

abr_IMB_CS = ncread([folder_path 'dIMB_CS.nc'], var_names);

abr_dLW_ts = ncread([folder_path 'dLW_ts.nc'], var_names);
abr_dLW_ta = ncread([folder_path 'dLW_ta.nc'], var_names);

abr_t = double(ncread([folder_path 'dGMST.nc'], 'year'));

%% Compute cloud feedback plus forcing masking

abr_dCLOUD_PLUS_FORCING = (abr_dLW_CS_PL - abr_dLW_PL) + (abr_dLW_CS_LR - abr_dLW_LR) + ...
    (abr_dSW_CS_AL - abr_dSW_AL) + (abr_dLW_CS_WV - abr_dLW_WV) + (abr_dSW_CS_WV - abr_dSW_WV) + ...
    (abr_IMB - abr_IMB_CS);

%% Compute radiative response and cloud feedback (from the fit!)

N = 3;
abr_beta_R = x(2*N+1:3*N);
abr_F_CS = x(9*N+1);
abr_F = sum(abr_beta_R);

abr_dCLOUD = abr_dCLOUD_PLUS_FORCING - (abr_F - abr_F_CS);
abr_dR = - (abr_IMB - abr_F);

%% Put in data

abr_Y_DATA = [abr_GMST'; abr_dR'; abr_dLW_PL'; abr_dLW_LR'; abr_dSW_AL'; abr_dLW_WV'; abr_dSW_WV'; abr_dCLOUD'];

%% Compute numerical green function

G_F = diff([ zeros(8,1), abr_Y_DATA], [], 2);

%% Compute radiative resposne
gamma = log(1.01)/log(4);
g_grad = gamma * t';
%g_grad = min(g_grad, ones(1,length(t'))); % Cap at 4xCO2, i.e. at gamma * t = 1

for j=1:length(t)
    resp_grad(:, j) = sum( G_F(:,1:j) .* g_grad(j:-1:1), 2 );
end

%% Append the plots


subplot(331)
plot(t, resp_grad(1,:), 'g-')

subplot(332)
plot(t, abr_F * gamma * t' - resp_grad(2,:), 'g-')

subplot(334)
plot(t, resp_grad(3,:), 'g-')

subplot(335)
plot(t, resp_grad(4,:), 'g-')

subplot(336)
plot(t, resp_grad(5,:), 'g-')

subplot(337)
plot(t, resp_grad(6,:), 'g-')

subplot(338)
plot(t, resp_grad(7,:), 'g-')

subplot(339)
plot(t, resp_grad(8,:), 'g-')

