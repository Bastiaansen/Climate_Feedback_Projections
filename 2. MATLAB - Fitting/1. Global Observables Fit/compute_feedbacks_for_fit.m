%% Extract parameters from fit
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

%% Compute contributions of different feedbacks
lambdas = - beta_R ./ beta_T;
lambdas_PL = beta_PL ./ beta_T;
lambdas_LR = beta_LR ./ beta_T;
lambdas_AL = beta_AL ./ beta_T;
lambdas_LWWV = beta_LWWV ./ beta_T;
lambdas_SWWV = beta_SWWV ./ beta_T;
lambdas_CLOUD = beta_CLOUD ./ beta_T;
lambdas_sum = lambdas_PL + lambdas_LR + lambdas_AL + lambdas_LWWV + lambdas_SWWV + lambdas_CLOUD;

%% Uncertainty propogation

% First derive the unidirectional uncertainty, which should be half of the
% total size of the confidence interval because of the assumed symmetry of
% the error distribution
uncertainty = (con(:,2)-con(:,1))/2;

% just the plain uncertainty
UN_tau = uncertainty(1:N);
UN_beta_T = uncertainty(N+1:2*N);
UN_beta_R = uncertainty(2*N+1:3*N);
UN_beta_PL = uncertainty(3*N+1:4*N);
UN_beta_LR = uncertainty(4*N+1:5*N);
UN_beta_AL = uncertainty(5*N+1:6*N);
UN_beta_LWWV = uncertainty(6*N+1:7*N);
UN_beta_SWWV = uncertainty(7*N+1:8*N);
UN_beta_CLOUD = uncertainty(8*N+1:9*N);
UN_F_CS = uncertainty(9*N+1);

% Propogation of uncertainy
% Computed via d(a/b)/|(a/b)| = da/|a| + db/|b|
% So d(a/b) = |a/b| * (da/|a| + db/|b|)

lambdas = - beta_R ./ beta_T;
lambdas_PL = beta_PL ./ beta_T;
lambdas_LR = beta_LR ./ beta_T;
lambdas_AL = beta_AL ./ beta_T;
lambdas_LWWV = beta_LWWV ./ beta_T;
lambdas_SWWV = beta_SWWV ./ beta_T;
lambdas_CLOUD = beta_CLOUD ./ beta_T;
lambdas_sum = lambdas_PL + lambdas_LR + lambdas_AL + lambdas_LWWV + lambdas_SWWV + lambdas_CLOUD;

UN_lambdas = propogate_confidenceinterval_quotient(lambdas,-beta_R,UN_beta_R,beta_T,UN_beta_T);
UN_lambdas_PL = propogate_confidenceinterval_quotient(lambdas,beta_PL,UN_beta_PL,beta_T,UN_beta_T);
UN_lambdas_LR = propogate_confidenceinterval_quotient(lambdas,beta_LR,UN_beta_LR,beta_T,UN_beta_T);
UN_lambdas_AL = propogate_confidenceinterval_quotient(lambdas,beta_AL,UN_beta_AL,beta_T,UN_beta_T);
UN_lambdas_LWWV = propogate_confidenceinterval_quotient(lambdas,beta_LWWV,UN_beta_LWWV,beta_T,UN_beta_T);
UN_lambdas_SWWV = propogate_confidenceinterval_quotient(lambdas,beta_SWWV,UN_beta_SWWV,beta_T,UN_beta_T);
UN_lambdas_CLOUD = propogate_confidenceinterval_quotient(lambdas,beta_CLOUD,UN_beta_CLOUD,beta_T,UN_beta_T);

% summing of uncertainties involves (assuming normal distributions) follows
% x sigma_f = \sqrt{ \sum (x \sigma_j)^2 }
UN_lambdas_sum = propogate_confidenceinterval_sum( [UN_lambdas_PL,...
    UN_lambdas_LR, UN_lambdas_AL, UN_lambdas_LWWV, UN_lambdas_SWWV, ...
    UN_lambdas_CLOUD ] );

% Propogation for sums and differences is the same.
UN_lambdas_mismatch = propogate_confidenceinterval_sum( [UN_lambdas, ...
    UN_lambdas_sum]);


%% display a report of the found feedback parameters.

%'±'

disp('------WHEN FITTING ' + string(length(tau)) + ' MODES ------------')
disp(' ')

for n = 1:length(lambdas)
    disp("MODE " + string(n) + ":")
    disp("tau (years) : " + string(tau(n)) + " ± " + string(UN_tau(n)))
    disp(' ')
    disp('lambda_n : ' + string(lambdas(n))+ " ± " + string(UN_lambdas(n)))
    disp(' ')
    disp('PLANCK fb: ' + string(lambdas_PL(n))+ " ± " + string(UN_lambdas_PL(n)))
    disp('LapseR fb: ' + string(lambdas_LR(n))+ " ± " + string(UN_lambdas_LR(n)))
    disp('ALBEDO fb: ' + string(lambdas_AL(n))+ " ± " + string(UN_lambdas_AL(n)))
    disp('WV LW  fb: ' + string(lambdas_LWWV(n))+ " ± " + string(UN_lambdas_LWWV(n)))
    disp('WV SW  fb: ' + string(lambdas_SWWV(n))+ " ± " + string(UN_lambdas_SWWV(n)))
    disp('CLOUD  fb: ' + string(lambdas_CLOUD(n))+ " ± " + string(UN_lambdas_CLOUD(n)))
    disp('sum fbs : ' + string(lambdas_sum(n))+ " ± " + string(UN_lambdas_sum(n)))
    disp('')
    disp('Mismatch: ' + string(lambdas_sum(n) - lambdas(n))+ " ± " + string(UN_lambdas_mismatch(n)))
    disp(' ')
    disp(' ')
end


%% Equilibrium feedback changes

R0 = sum(beta_R);
T_EQ = sum(beta_T);
PL_EQ = sum(beta_PL);
LR_EQ = sum(beta_LR);
AL_EQ = sum(beta_AL);
LWWV_EQ = sum(beta_LWWV);
SWWV_EQ = sum(beta_SWWV);
CLOUD_EQ = sum(beta_CLOUD);

lambda_EQ = - R0 / T_EQ;

lambda_PL_EQ = PL_EQ / T_EQ;
lambda_LR_EQ = LR_EQ / T_EQ;
lambda_AL_EQ = AL_EQ / T_EQ;
lambda_LWWV_EQ = LWWV_EQ / T_EQ;
lambda_SWWV_EQ = SWWV_EQ / T_EQ;
lambda_CLOUD_EQ = CLOUD_EQ / T_EQ;

lambda_sum_EQ = lambda_PL_EQ + lambda_LR_EQ + lambda_AL_EQ + lambda_LWWV_EQ + lambda_SWWV_EQ + lambda_CLOUD_EQ;

%% Equilibrium uncertainties (i.e. propogate uncertainties)

UN_R0 =  propogate_confidenceinterval_sum( UN_beta_R' );
UN_T_EQ = propogate_confidenceinterval_sum( UN_beta_T' );
UN_PL_EQ = propogate_confidenceinterval_sum( UN_beta_PL' );
UN_LR_EQ = propogate_confidenceinterval_sum( UN_beta_LR' );
UN_AL_EQ = propogate_confidenceinterval_sum( UN_beta_AL' );
UN_LWWV_EQ = propogate_confidenceinterval_sum( UN_beta_LWWV' );
UN_SWWV_EQ = propogate_confidenceinterval_sum( UN_beta_SWWV' );
UN_CLOUD_EQ = propogate_confidenceinterval_sum( UN_beta_CLOUD' );

UN_lambda_EQ =  propogate_confidenceinterval_quotient(lambda_EQ,-R0,UN_R0,T_EQ,UN_T_EQ);

UN_lambda_PL_EQ = propogate_confidenceinterval_quotient(lambda_PL_EQ,PL_EQ,UN_PL_EQ,T_EQ,UN_T_EQ);
UN_lambda_LR_EQ = propogate_confidenceinterval_quotient(lambda_LR_EQ,LR_EQ,UN_LR_EQ,T_EQ,UN_T_EQ);
UN_lambda_AL_EQ = propogate_confidenceinterval_quotient(lambda_AL_EQ,AL_EQ,UN_AL_EQ,T_EQ,UN_T_EQ);
UN_lambda_LWWV_EQ = propogate_confidenceinterval_quotient(lambda_LWWV_EQ,LWWV_EQ,UN_LWWV_EQ,T_EQ,UN_T_EQ);
UN_lambda_SWWV_EQ = propogate_confidenceinterval_quotient(lambda_SWWV_EQ,SWWV_EQ,UN_SWWV_EQ,T_EQ,UN_T_EQ);
UN_lambda_CLOUD_EQ = propogate_confidenceinterval_quotient(lambda_CLOUD_EQ,CLOUD_EQ,UN_CLOUD_EQ,T_EQ,UN_T_EQ);

UN_lambda_sum_EQ = propogate_confidenceinterval_sum( [ UN_lambda_PL_EQ, ...
    UN_lambda_LR_EQ, UN_lambda_AL_EQ, UN_lambda_LWWV_EQ, UN_lambda_SWWV_EQ, ...
    UN_lambda_CLOUD_EQ ] );

UN_lambda_mismatch_EQ = propogate_confidenceinterval_sum( [ UN_lambda_EQ, ...
    UN_lambda_sum_EQ ] );

%% display a report of the found feedback parameters in equilibrium


disp("EQUILIBRIUM")
disp("tau (years) : N/A")
disp(' ')
disp("lambda_EQ : " + string(lambda_EQ) + " ± " + string(UN_lambda_EQ))
disp(' ')
disp('PLANCK fb: ' + string(lambda_PL_EQ) + " ± " + string(UN_lambda_PL_EQ))
disp('LapseR fb: ' + string(lambda_LR_EQ) + " ± " + string(UN_lambda_LR_EQ))
disp('ALBEDO fb: ' + string(lambda_AL_EQ) + " ± " + string(UN_lambda_AL_EQ))
disp('WV LW  fb: ' + string(lambda_LWWV_EQ) + " ± " + string(UN_lambda_LWWV_EQ))
disp('WV SW  fb: ' + string(lambda_SWWV_EQ) + " ± " + string(UN_lambda_SWWV_EQ))
disp('CLOUD  fb: ' + string(lambda_CLOUD_EQ) + " ± " + string(UN_lambda_CLOUD_EQ))
disp('sum fbs : ' + string(lambda_sum_EQ) + " ± " + string(UN_lambda_sum_EQ))
disp('')
disp('Mismatch: ' + string(lambda_sum_EQ - lambda_EQ) + " ± " + string(UN_lambda_mismatch_EQ))
disp(' ')
disp(' ')