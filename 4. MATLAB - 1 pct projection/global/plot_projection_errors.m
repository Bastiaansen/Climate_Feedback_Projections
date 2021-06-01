
figure('Units','normalized','Position', [0.05 0 0.9 0.9]);


Y_DATA = [GMST'; IMB'; dLW_PL'; dLW_LR'; dSW_AL'; dLW_WV'; dSW_WV'; dCLOUD_PLUS_FORCING'];
error_Y = (Y(:,1:150) - Y_DATA);

subplot(331)
hold on
plot(t, error_Y(1,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('GMST')

subplot(332)
hold on
plot(t, error_Y(2,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('IMB')

subplot(334)
hold on
plot(t, error_Y(3,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dLW PL')

subplot(335)
hold on
plot(t, error_Y(4,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dLW LR')

subplot(336)
hold on
plot(t, error_Y(5,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dSW AL')

subplot(337)
hold on
plot(t, error_Y(6,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dLW WV')

subplot(338)
hold on
plot(t, error_Y(7,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dSW WV')

subplot(339)
hold on
plot(t, error_Y(8,:), 'Color', 'k')
xlabel('$t$', 'Interpreter', 'latex')
title('dCLOUD plus Forcing')
