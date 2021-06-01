% This script automatically finds the computed observable files if left in
% the original folders. Then it computes cloud feedback contributions and
% tries varies values of amount of modes (M in the paper, but N in this
% script) from M = 1 to M = 5 and creates some statistics from these tests.
% It also creates plots of best fits for different M values. Finally, M = 3
% is taken to create the plots in the paper and to compute feedback
% strengths per mode.

%% clean slate

clear all
close all

%% Load-in data
var_names = '__xarray_dataarray_variable__';
folder_path = '..\..\1. Python - CMIP Feedback Computations\CESM2 abrupt4xCO2\Data\global\';

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



%% Try different values for the amount of modes and initial guesses in the regression

dCLOUD_PLUS_FORCING = (dLW_CS_PL - dLW_PL) + (dLW_CS_LR - dLW_LR) + ...
    (dSW_CS_AL - dSW_AL) + (dLW_CS_WV - dLW_WV) + (dSW_CS_WV - dSW_WV) + ...
    (IMB - IMB_CS);

Y_DATA = [GMST'; IMB'; dLW_PL'; dLW_LR'; dSW_AL'; dLW_WV'; dSW_WV'; dCLOUD_PLUS_FORCING'];


maxN = 5; % The maximum number of modes to consider
M = 1000; % How many initial guesses to test

x_fits_all = cell(maxN,M);
ress = zeros(maxN,M);
confidence_intervals_all = cell(maxN,M);

for i=1:maxN
    tau = linspace(1,1000,i)';    
    N = length(tau)

    for j=1:M
    
        beta_T = 5*randn(N,1);
        beta_R = 5*randn(N,1);
        beta_PL = -5*randn(N,1);
        beta_LR = -5*randn(N,1);
        beta_AL = 5*randn(N,1);
        beta_LWWV = 5*randn(N,1);
        beta_SWWV = 5*randn(N,1);
        beta_CLOUD = 5*randn(N,1);
        F_CS = randn(1,1);

        x_init = [tau; beta_T; beta_R; beta_PL; beta_LR; beta_AL; beta_LWWV; beta_SWWV; beta_CLOUD; F_CS];

        opts1 = optimset('display', 'off');
        [x_fit, resnorm, residual, exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,t) exponent_fit_function(x,t), x_init, t', Y_DATA, [], [], opts1);
 
        conf = nlparci(x_fit,residual,'jacobian',jacobian);
                
        ress(i,j) = resnorm;
        x_fits_all{i,j} = x_fit;
        confidence_intervals_all{i,j} = conf;
        
    end

end


%% Statistics on fits with different initial conditions
% Compute mean values, 5 and 95 percentile values as well as best values
% (all measured via the total residual).

ress_mean = mean(ress,2);
ress_5 =  prctile(ress,5,2);
ress_95 =  prctile(ress,95,2);
[ress_best, ress_best_index] = min(ress,[],2);

% for plotting purposes we will use the best fits only
x_fits = {};
confidence_intervals = {};

for i=1:maxN
    x_fits{i} = x_fits_all{i,ress_best_index(i)};
    confidence_intervals{i} = confidence_intervals_all{i,ress_best_index(i)};
end

%% Save data
fit_data_folder = 'fit_data';
% Save all fits:
save([fit_data_folder '/all_fits.mat'], 'ress', 'x_fits_all', 'confidence_intervals_all');
% Save the best fits
save([fit_data_folder '/best_fits.mat'], 'x_fits', 'confidence_intervals')

%% Load data
fit_data_folder = 'fit_data';
load([fit_data_folder '/all_fits.mat'])
load([fit_data_folder '/best_fits.mat'])

%% Plotting + FEEDBACK COMPUTATION


figure('Units','normalized','Position', [0.05 0 0.9 0.9]);

subplot(331)
plot(t, GMST, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('GMST')

subplot(332)
plot(t,IMB, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('IMB')

subplot(334)
plot(t,dLW_PL, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW PL')

subplot(335)
plot(t, dLW_LR, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW LR')

subplot(336)
plot(t, dSW_AL, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dSW AL')

subplot(337)
plot(t, dLW_WV, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW WV')

subplot(338)
plot(t, dSW_WV, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dSW WV')

subplot(339)
plot(t, dCLOUD_PLUS_FORCING, 'ro', 'MarkerSize', 1, 'HandleVisibility', 'off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dCLOUD plus Forcing')


colors = jet(length(x_fits));

for i = 1:length(x_fits)
    
    x = x_fits{i};
    
    append_plot_with_fit
    
end

legend(string([1:maxN]))


subplot(332)
ylim([0 10])


%% Plots of different regression methods and their residuals

figure('Units','normalized','Position', [0.05 0 0.9 0.9]);

bar_data = [ ress_5, ress_mean, ress_95];

bar(1:maxN, bar_data, 1)

hold on
scatter(1:maxN, ress_best, 200, 'k*')
xlabel('$M$', 'Interpreter', 'latex', 'fontsize', 20)
ylabel('Residual', 'Interpreter', 'latex', 'fontsize', 20)

%% Display the fitted/found timescales for all N
for i=1:length(x_fits)
    x_fiti = x_fits{i};
    N = (length(x_fiti)-1)/9;
    tauI = x_fiti(1:N)
end

%% Computation of feedback for N = 3 (i.e. results in paper)
x = x_fits{3};
con = confidence_intervals{3};
compute_feedbacks_for_fit

%% Figure with only plot for N = 3
figure('Units','normalized','Position', [0.05 0 0.9 0.9]);

subplot(331)
plot(t, GMST, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('GMST')

subplot(332)
plot(t,IMB, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('IMB')

subplot(334)
plot(t,dLW_PL, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW PL')

subplot(335)
plot(t, dLW_LR, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW LR')

subplot(336)
plot(t, dSW_AL, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dSW AL')

subplot(337)
plot(t, dLW_WV, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dLW WV')

subplot(338)
plot(t, dSW_WV, 'ro', 'MarkerSize', 1, 'HandleVisibility','off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dSW WV')

F_CS = x(9*N+1);
F_0 = sum(beta_R);

subplot(339)
plot(t, dCLOUD_PLUS_FORCING - (F_0 - F_CS), 'ro', 'MarkerSize', 1, 'HandleVisibility', 'off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dCLOUD')

i = 1;
colors = [0.1 0.1 1];

% Do some 'magic' with variables to just use t for the fits, but we want to
% make those fits for more than 1,000 years so we need to use a different
% variable t. At the end revert back to the time points for the data.
t_data = t;
t_fit = [0:1:3000]';
t = t_fit;
append_plot_with_fit_CLOUD
t = t_data;
