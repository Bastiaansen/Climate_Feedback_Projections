% This script plots the data for the 1pct scenario along side the
% projections made à la Linear Response Theory with a data Green Function
% and with a Green Function obtained from nonlinear fits.

%% First: Obtain estimations for forcing and cloud masking of forcing
beta_factor = log(1.01) / log(4);
N = (length(x)-1)/9;
F_FS = sum(x(2*N+1:3*N)) .* beta_factor;
F_CS = x(9*N+1) .* beta_factor;

F_cloud_masking_DATA = (F_FS-F_CS) * t';
F_cloud_masking_fit = (F_FS-F_CS) * t_fit';

%%
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
plot(t, dCLOUD_PLUS_FORCING - F_cloud_masking_DATA, 'ro', 'MarkerSize', 1, 'HandleVisibility', 'off')
hold on
xlabel('$t$', 'Interpreter', 'latex')
title('dCLOUD')



%% Add the projections


subplot(331)
plot(t_fit, Y(1,:), 'Color', 'b')

subplot(332)
plot(t_fit,Y(2,:), 'Color', 'b')

subplot(334)
plot(t_fit,Y(3,:), 'Color', 'b')

subplot(335)
plot(t_fit, Y(4,:), 'Color', 'b')

subplot(336)
plot(t_fit, Y(5,:), 'Color', 'b')

subplot(337)
plot(t_fit, Y(6,:), 'Color', 'b')

subplot(338)
plot(t_fit, Y(7,:), 'Color', 'b')

subplot(339)
plot(t_fit, Y(8,:) - F_cloud_masking_fit, 'Color', 'b')

%% Add projections from linear response theory

add_LRT_projections