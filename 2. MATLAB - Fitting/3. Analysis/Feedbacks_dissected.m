% This script should be run while in the folder with the exponent_fit_function,
% i.e. the folder '1. Global Observables Fit', as it relies on some
% functions located in that folder.

% This script uses the found fitted parameters to compute the instantaneous
% feedback change (as the local slope of (DT,DRj).

%% Load-in fitted parameters
fit_data_folder = '../../2. MATLAB - Fitting/1. Global Observables Fit/fit_data';
load([fit_data_folder '/best_fits.mat'])
% We use the M = 3 scenario here
M = 3;
x = x_fits{3};

%% Set-up for plots and compute the dynamics for all observables


T_end = 10^5;
T_end = 1000;
t = 0:1:T_end;
Y = exponent_fit_function(x, t);

%% Extract the different observables
N = (length(x)-1)/9;
F_0 = sum(x(2*N+1:3*N));
F_CS = x(9*N+1);

DT = Y(1,:);
N = Y(2,:);
R = N - F_0;
R_PL = Y(3,:);
R_LR = Y(4,:);
R_SA = Y(5,:);
R_WV_LW = Y(6,:);
R_WV_SW = Y(7,:);
R_CL = Y(8,:) - (F_0 - F_CS);

%% Plot evolution of observables on semilog scale

figure('Units','normalized','Position', [0.05 0 0.9 0.9]);
semilogx(t, R, 'k--')
hold on
semilogx(t, R_PL)
semilogx(t, R_LR)
semilogx(t, R_SA)
semilogx(t, R_WV_LW)
semilogx(t, R_WV_SW)
semilogx(t, R_CL)
xlim([0 T_end])
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$R_j$', 'Interpreter', 'latex')
legend({'$R$', '$R_{PL}$', '$R_{LR}$', '$R_{SA}$', '$R_{WV LW}$', '$R_{WV SW}$', '$R_{CL}$'}, 'Interpreter', 'latex')

%% Plot the fraction/contribution per observable on semilog scale
figure('Units','normalized','Position', [0.05 0 0.9 0.9]);
hold on
semilogx(t, R./R, 'k--')
hold on
semilogx(t, R_PL./R)
semilogx(t, R_LR./R)
semilogx(t, R_SA./R)
semilogx(t, R_WV_LW./R)
semilogx(t, R_WV_SW./R)
semilogx(t, R_CL./R)
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$R_j$', 'Interpreter', 'latex')
legend({'$R$', '$R_{PL}$', '$R_{LR}$', '$R_{SA}$', '$R_{WV LW}$', '$R_{WV SW}$', '$R_{CL}$'}, 'Interpreter', 'latex')
%% Plot the incremental change in the feedback contribution for the abrupt4xCO2 scenario

Y_der = exponent_fit_function_der(x, t);

DT_der = Y_der(1,:);
DR_der = Y_der(2,:);
DR_PL_der = Y_der(3,:);
DR_LR_der = Y_der(4,:);
DR_SA_der = Y_der(5,:);
DR_WV_LW_der = Y_der(6,:);
DR_WV_SW_der = Y_der(7,:);
DR_CL = Y_der(8,:);

lambda_R2 = DR_der ./ DT_der;
lambda_R2_PL = DR_PL_der ./ DT_der;
lambda_R2_LR = DR_LR_der ./ DT_der;
lambda_R2_SA = DR_SA_der ./ DT_der;
lambda_R2_WV_LW = DR_WV_LW_der ./ DT_der;
lambda_R2_WV_SW = DR_WV_SW_der ./ DT_der;
lambda_R2_CL = DR_CL ./ DT_der;

lambda_R2_mismatch = (lambda_R2_PL + lambda_R2_LR + lambda_R2_SA + lambda_R2_WV_LW + lambda_R2_WV_SW + lambda_R2_CL) - lambda_R2;


figure('Units','normalized','Position', [0.05 0 0.9 0.9]);
plot(t, lambda_R2, 'k:')
hold on
semilogx(t, lambda_R2_PL)
semilogx(t, lambda_R2_LR)
semilogx(t, lambda_R2_SA)
semilogx(t, lambda_R2_WV_LW)
semilogx(t, lambda_R2_WV_SW)
semilogx(t, lambda_R2_CL)
semilogx(t, lambda_R2_mismatch, 'r--')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\tilde\lambda_j$', 'Interpreter', 'latex')
legend({'$\lambda$', '$\lambda_{PL}$', '$\lambda_{LR}$', '$\lambda_{SA}$', '$\lambda_{WV LW}$', '$\lambda_{WV SW}$', '$\lambda_{CL}$', 'remaining'}, 'Interpreter', 'latex')
title('Feedback strengths (local slope) for abupt4xCO2')
%% Plot the incremental change in the feedback contribution for the 1pctCO2 scenario
% We do not need to compute any additional things, since the derivatives of
% the observable evolution in the 1pctCO2 scenario correspond with the
% evolution of the observations in the abrupt4xCO2 scenario. Note that even
% the factor log(1.01)/log(4) that you need to relate 1pctCO2 projections
% to the abrupt4xCO2 fits is not necessary here, since we look at a
% fraction of things that both have this factor.

lambda_R = R ./ DT;
lambda_R_PL = R_PL ./ DT;
lambda_R_LR = R_LR ./ DT;
lambda_R_SA = R_SA ./ DT;
lambda_R_WV_LW = R_WV_LW ./ DT;
lambda_R_WV_SW = R_WV_SW ./ DT;
lambda_R_CL = R_CL ./ DT;

lambda_R_mismatch = (lambda_R_PL + lambda_R_LR + lambda_R_SA + lambda_R_WV_LW + lambda_R_WV_SW + lambda_R_CL) - lambda_R;

figure('Units','normalized','Position', [0.05 0 0.9 0.9]);
plot(t, lambda_R, 'k:')
hold on
semilogx(t, lambda_R_PL)
semilogx(t, lambda_R_LR)
semilogx(t, lambda_R_SA)
semilogx(t, lambda_R_WV_LW)
semilogx(t, lambda_R_WV_SW)
semilogx(t, lambda_R_CL)
semilogx(t, lambda_R_mismatch, 'r--')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\tilde\lambda_j$', 'Interpreter', 'latex')
legend({'$\lambda$', '$\lambda_{PL}$', '$\lambda_{LR}$', '$\lambda_{SA}$', '$\lambda_{WV LW}$', '$\lambda_{WV SW}$', '$\lambda_{CL}$', 'remaining'}, 'Interpreter', 'latex')
title('Feedback strengths (local slope) for 1pctCO2')

