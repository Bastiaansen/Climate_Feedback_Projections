%% EXTRACT MODEL PARAMETERS

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


dT_eq = sum(beta_T);
F_0 = sum(beta_R);
dPL_eq = sum(beta_PL);
dLR_eq = sum(beta_LR);
dAL_eq = sum(beta_AL);
dLWWV_eq = sum(beta_LWWV);
dSWWV_eq = sum(beta_SWWV);
dCLOUD_eq = sum(beta_CLOUD);

Y = exponent_fit_function(x,t');

%% PLOTTING

col = colors(i,:);

subplot(331)
plot(t, Y(1,:), 'Color', col)
plot([0 t(end)], dT_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(332)
plot(t,Y(2,:), 'Color', col)
plot([0 t(end)], F_0 * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(334)
plot(t,Y(3,:), 'Color', col)
plot([0 t(end)], dPL_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(335)
plot(t, Y(4,:), 'Color', col)
plot([0 t(end)], dLR_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(336)
plot(t, Y(5,:), 'Color', col)
plot([0 t(end)], dAL_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(337)
plot(t, Y(6,:), 'Color', col)
plot([0 t(end)], dLWWV_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(338)
plot(t, Y(7,:), 'Color', col)
plot([0 t(end)], dSWWV_eq * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')

subplot(339)
plot(t, Y(8,:), 'Color', col)
plot([0 t(end)], dCLOUD_eq + (F_0 - F_CS) * [1,1], 'Color', col, 'LineStyle', '--', 'HandleVisibility','off')