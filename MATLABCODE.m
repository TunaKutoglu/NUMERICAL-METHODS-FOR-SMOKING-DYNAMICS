% Parametreler
A = 0.25; beta = 0.03; d = 0.004; mu = 0.008;
epsilon = 0.02; delta = 0.074; alpha = 0.15; b = 0.4;
% Adım büyüklüğü ve simülasyon süresi
h = 0.1;
t_end = 500;
t = 0:h:t_end;
N = length(t);
% Başlangıç değer kombinasyonları
initial_conditions = [
    100, 30, 10, 20; % İlk başlangıç seti
    30, 10, 20, 100; % İkinci başlangıç seti
    20, 100, 30, 10  % Üçüncü başlangıç seti
];
% Her başlangıç seti için grafikler
for i = 1:size(initial_conditions, 1)
    % Başlangıç değerlerini al
    P0 = initial_conditions(i, 1);
    L0 = initial_conditions(i, 2);
    S0 = initial_conditions(i, 3);
    Q0 = initial_conditions(i, 4);
    % Euler, Modified Euler ve Runge-Kutta için başlangıç koşulları
    P_euler = zeros(1, N); L_euler = zeros(1, N); S_euler = zeros(1, N); Q_euler = zeros(1, N);
    P_mod = zeros(1, N); L_mod = zeros(1, N); S_mod = zeros(1, N); Q_mod = zeros(1, N);
    P_rk2 = zeros(1, N); L_rk2 = zeros(1, N); S_rk2 = zeros(1, N); Q_rk2 = zeros(1, N);
    P_euler(1) = P0; L_euler(1) = L0; S_euler(1) = S0; Q_euler(1) = Q0;
    P_mod(1) = P0; L_mod(1) = L0; S_mod(1) = S0; Q_mod(1) = Q0;
    P_rk2(1) = P0; L_rk2(1) = L0; S_rk2(1) = S0; Q_rk2(1) = Q0;
    % Sayısal yöntemlerin uygulanması
    for n = 1:N-1
        % Euler yöntemi
        dP_euler = A - (2 * beta * P_euler(n) * L_euler(n)) / (P_euler(n) + L_euler(n)) - mu * P_euler(n);
        dL_euler = (2 * beta * P_euler(n) * L_euler(n)) / (P_euler(n) + L_euler(n)) - (epsilon + d + mu) * L_euler(n) + b * alpha * Q_euler(n);
        dS_euler = epsilon * L_euler(n) - delta * S_euler(n) - (d + mu) * S_euler(n) - (1 - b) * alpha * Q_euler(n);
        dQ_euler = delta * S_euler(n) - (mu + d) * Q_euler(n) - alpha * Q_euler(n);
        P_euler(n+1) = P_euler(n) + h * dP_euler;
        L_euler(n+1) = L_euler(n) + h * dL_euler;
        S_euler(n+1) = S_euler(n) + h * dS_euler;
        Q_euler(n+1) = Q_euler(n) + h * dQ_euler;
        % Modified Euler yöntemi
        dP_mod = dP_euler; dL_mod = dL_euler; dS_mod = dS_euler; dQ_mod = dQ_euler;
        P_star = P_mod(n) + h * dP_mod;
        L_star = L_mod(n) + h * dL_mod;
        S_star = S_mod(n) + h * dS_mod;
        Q_star = Q_mod(n) + h * dQ_mod;
        dP_star = A - (2 * beta * P_star * L_star) / (P_star + L_star) - mu * P_star;
        dL_star = (2 * beta * P_star * L_star) / (P_star + L_star) - (epsilon + d + mu) * L_star + b * alpha * Q_star;
        dS_star = epsilon * L_star - delta * S_star - (d + mu) * S_star - (1 - b) * alpha * Q_star;
        dQ_star = delta * S_star - (mu + d) * Q_star - alpha * Q_star;
        P_mod(n+1) = P_mod(n) + (h / 2) * (dP_mod + dP_star);
        L_mod(n+1) = L_mod(n) + (h / 2) * (dL_mod + dL_star);
        S_mod(n+1) = S_mod(n) + (h / 2) * (dS_mod + dS_star);
        Q_mod(n+1) = Q_mod(n) + (h / 2) * (dQ_mod + dQ_star);
        % 2nd Order Runge-Kutta yöntemi
        k1_P = dP_euler; k1_L = dL_euler; k1_S = dS_euler; k1_Q = dQ_euler;
        P_star = P_rk2(n) + h * k1_P;
        L_star = L_rk2(n) + h * k1_L;
        S_star = S_rk2(n) + h * k1_S;
        Q_star = Q_rk2(n) + h * k1_Q;
        k2_P = A - (2 * beta * P_star * L_star) / (P_star + L_star) - mu * P_star;
        k2_L = (2 * beta * P_star * L_star) / (P_star + L_star) - (epsilon + d + mu) * L_star + b * alpha * Q_star;
        k2_S = epsilon * L_star - delta * S_star - (d + mu) * S_star - (1 - b) * alpha * Q_star;
        k2_Q = delta * S_star - (mu + d) * Q_star - alpha * Q_star;
        P_rk2(n+1) = P_rk2(n) + (h / 2) * (k1_P + k2_P);
        L_rk2(n+1) = L_rk2(n) + (h / 2) * (k1_L + k2_L);
        S_rk2(n+1) = S_rk2(n) + (h / 2) * (k1_S + k2_S);
        Q_rk2(n+1) = Q_rk2(n) + (h / 2) * (k1_Q + k2_Q);
    end
    % Grafikler
    figure;
    subplot(2, 2, 1);
    plot(t, P_euler, 'b--', t, P_mod, 'r-.', t, P_rk2, 'g:', 'LineWidth', 1.5);
    xlabel('t'); ylabel('P(t)'); title('P(t)');
    legend('Euler', 'Modified Euler', '2nd Order Runge-Kutta'); grid on;
    subplot(2, 2, 2);
    plot(t, L_euler, 'b--', t, L_mod, 'r-.', t, L_rk2, 'g:', 'LineWidth', 1.5);
    xlabel('t'); ylabel('L(t)'); title('L(t)');
    legend('Euler', 'Modified Euler', '2nd Order Runge-Kutta'); grid on;
    subplot(2, 2, 3);
    plot(t, S_euler, 'b--', t, S_mod, 'r-.', t, S_rk2, 'g:', 'LineWidth', 1.5);
    xlabel('t'); ylabel('S(t)'); title('S(t)');
    legend('Euler', 'Modified Euler', '2nd Order Runge-Kutta'); grid on;
    subplot(2, 2, 4);
    plot(t, Q_euler, 'b--', t, Q_mod, 'r-.', t, Q_rk2, 'g:', 'LineWidth', 1.5);
    xlabel('t'); ylabel('Q(t)'); title('Q(t)');
    legend('Euler', 'Modified Euler', '2nd Order Runge-Kutta'); grid on;
    sgtitle(['P0=', num2str(P0), ', L0=', num2str(L0), ', S0=', num2str(S0), ', Q0=', num2str(Q0)]);
end
