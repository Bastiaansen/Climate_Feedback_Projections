function Y = exponent_fit_function(x,t)

N = (length(x)-1)/9;

tau = x(1:N);
beta_T = x(N+1:2*N);
beta_R = x(2*N+1:3*N);
beta_PL = x(3*N+1:4*N);
beta_LR = x(4*N+1:5*N);
beta_AL = x(5*N+1:6*N);
beta_LWWV = x(6*N+1:7*N);
beta_SWWV = x(7*N+1:8*N);
beta_CLOUD = x(8*N+1:9*N);
F_CS = x(9*N+1);

Modes = 1 - exp(- repmat(t,length(tau),1) ./ tau );

dT = beta_T' * Modes;
dR = sum(beta_R) - beta_R' * Modes;
dPL = beta_PL' * Modes;
dLR = beta_LR' * Modes;
dAL = beta_AL' * Modes;
dLWWV = beta_LWWV' * Modes;
dSWWV = beta_SWWV' * Modes;
dCLOUD_PLUS_FORCING = ( sum(beta_R) - F_CS ) + beta_CLOUD' * Modes;


Y = [dT; dR; dPL; dLR; dAL; dLWWV; dSWWV; dCLOUD_PLUS_FORCING];


end

