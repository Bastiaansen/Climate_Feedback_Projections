function Y = exponent_fit_function_1pct(x,t)
% The fits have been made for the abrupt-4xCO2 experiment. To apply them to
% the 1pct experiment we have to transform  the coefficients, as in the
% 4xCO2 experiment we have tilde_beta_m = beta_m * A * log(4) and in the 
% 1pctCO2 experimetn we have hat_beta_m = beta_m * A * log(1.01). Hence
% transforming the tilde_beta values (that have been fitted) to the
% hat_beta_m values can be done by multiplying with log(1.01)/log(4), which
% is the beta_factor introduced below.
% In the 1pct situation, the modes are
% hat_beta_m * (t + tau_m * (exp(-t/tau_m)-1))

N = (length(x)-1)/9;

tau = x(1:N);

beta_factor = log(1.01) / log(4);

beta_T = x(N+1:2*N) .* beta_factor;
beta_R = x(2*N+1:3*N) .* beta_factor;
beta_PL = x(3*N+1:4*N) .* beta_factor;
beta_LR = x(4*N+1:5*N) .* beta_factor;
beta_AL = x(5*N+1:6*N) .* beta_factor;
beta_LWWV = x(6*N+1:7*N) .* beta_factor;
beta_SWWV = x(7*N+1:8*N) .* beta_factor;
beta_CLOUD = x(8*N+1:9*N) .* beta_factor;
F_CS = x(9*N+1);

Modes = repmat(t,length(tau),1) + tau .* ( exp(- repmat(t,length(tau),1) ./ abs(tau)) - 1);

dT = beta_T' * Modes;
dR = sum(beta_R) * t - beta_R' * Modes;
dPL = beta_PL' * Modes;
dLR = beta_LR' * Modes;
dAL = beta_AL' * Modes;
dLWWV = beta_LWWV' * Modes;
dSWWV = beta_SWWV' * Modes;
dCLOUD_PLUS_FORCING = ( sum(beta_R) - F_CS * beta_factor) * t + beta_CLOUD' * Modes;


Y = [dT; dR; dPL; dLR; dAL; dLWWV; dSWWV; dCLOUD_PLUS_FORCING];

end

