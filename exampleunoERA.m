% =========================================================================
% MATLAB Script: ECO-Climatic plant--pollinator (Canonical Scenario A) 
% date creation: 5/oct/2026
% date modification: 25/nov/2026
% Author: Anibal Coronel
% =========================================================================

clear; 
clc; 
close all;

% --- Time interval for simulation (15-year horizon, S4-6) ---
t_span = [0 15]; 
omega  = 1.0;

% Solution 1: Transformed Low initial conditions (1/x0, 1/y0, C0, T0)
IC1 = [1/1.5, 1/1.0, 1.0, 12.0]; 

% Solution 2: Transformed High initial conditions (1/x0, 1/y0, C0, T0)
IC2 = [1/12.0, 1/10.0, 8.0, 20.0]; 

% --- Numerical solution with ODE45 ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Compute exact omega-periodic fixed point u*(0) via Newton shooting (S4-5)
u_star0 = [1/9.13; 1/3.96; 1.98; 16.5];
for iter = 1:15
    [~, U_per] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
        [0 omega], u_star0, options);
    F_shoot = U_per(end, :)' - u_star0;
    if norm(F_shoot, 2) < 1e-11
        break;
    end
    J_shoot = zeros(4, 4);
    eps_fd  = 1e-7;
    for j = 1:4
        du = zeros(4, 1); du(j) = eps_fd;
        [~, U_p] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
            [0 omega], u_star0 + du, options);
        J_shoot(:, j) = (U_p(end, :)' - U_per(end, :)') / eps_fd ...
            - (j == (1:4)');
    end
    u_star0 = u_star0 - J_shoot \ F_shoot;
end

[t1, sol1] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
    t_span, IC1, options);
[t2, sol2] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
    t_span, IC2, options);
[ts, sols] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
    t_span, u_star0, options);

% --- Figure 3 generation (preserving original fonts & line widths) ---
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

subplot(2, 2, 1);
plot(t1, sol1(:, 1), 'r-', t2, sol2(:, 1), 'b:', 'LineWidth', 3);
hold on;
plot(ts, sols(:, 1), 'k--', 'LineWidth', 2.0);
xlim([0 15]);
legend('Low i.c.', 'High i.c.', 'Exact $\mathbf{u}^*(t)$', ...
    'FontSize', 13, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_1(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Plants $u_1(t) = 1/x(t)$', ...
    'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 2);
plot(t1, sol1(:, 2), 'r-', t2, sol2(:, 2), 'b:', 'LineWidth', 3);
hold on;
plot(ts, sols(:, 2), 'k--', 'LineWidth', 2.0);
xlim([0 15]);
legend('Low i.c.', 'High i.c.', 'Exact $\mathbf{u}^*(t)$', ...
    'FontSize', 13, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_2(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Pollinators $u_2(t) = 1/y(t)$', ...
    'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 3);
plot(t1, sol1(:, 3), 'r-', t2, sol2(:, 3), 'b:', 'LineWidth', 3);
hold on;
plot(ts, sols(:, 3), 'k--', 'LineWidth', 2.0);
xlim([0 15]);
legend('Low i.c.', 'High i.c.', 'Exact $\mathbf{u}^*(t)$', ...
    'FontSize', 13, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_3(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Greenhouse Gases $u_3(t) = C(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 4);
plot(t1, sol1(:, 4), 'r-', t2, sol2(:, 4), 'b:', 'LineWidth', 3);
hold on;
plot(ts, sols(:, 4), 'k--', 'LineWidth', 2.0);
xlim([0 15]);
legend('Low i.c.', 'High i.c.', 'Exact $\mathbf{u}^*(t)$', ...
    'FontSize', 13, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_4(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Temperature Elevation $u_4(t) = T(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

% --- Zoom insets for u1, u2, u3 on [10, 15] x [0.05, 3] (S4-7) ---
idx1 = (t1 >= 10 & t1 <= 15);
idx2 = (t2 >= 10 & t2 <= 15);
idxs = (ts >= 10 & ts <= 15);

inset_pos = [
    0.22, 0.68, 0.18, 0.12;   % Subplot (1) u1
    0.66, 0.68, 0.18, 0.12;   % Subplot (2) u2
    0.22, 0.23, 0.18, 0.12    % Subplot (3) u3:
];

for k = 1:3
    ax_in = axes('Position', inset_pos(k, :));
    box(ax_in, 'on');
    plot(ax_in, t1(idx1), sol1(idx1, k), 'r-', ...
        t2(idx2), sol2(idx2, k), 'b:', 'LineWidth', 2.5);
    hold(ax_in, 'on');
    plot(ax_in, ts(idxs), sols(idxs, k), 'k--', 'LineWidth', 1.8);
    xlim(ax_in, [10 15]);
    ylim(ax_in, [0.05 3.0]);
    title(ax_in, 'Zoom: $t \in [10, 15]$', ...
        'FontSize', 11, 'Interpreter', 'latex');
    set(ax_in, 'FontName', 'Times New Roman', 'FontSize', ...
        11, 'FontWeight', 'bold', 'LineWidth', 1.5);
    grid(ax_in, 'on');
end

% --- Exportation to folder ---
folder_name = 'FIGURASERA';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

file_path_png = fullfile(folder_name, 'convergenia_periodic_orbit.png');
file_path_eps = fullfile(folder_name, 'convergenia_periodic_orbit.eps');

set(fig, 'PaperPositionMode', 'auto');
print(fig, file_path_png, '-dpng', '-r300'); 
print(fig, file_path_eps, '-depsc', '-r300'); 

% --- Periodicity Residual Diagnostic across ODE Tolerances (S4-6) ---
tols   = [1e-4, 1e-6, 1e-8, 1e-10];
styles = {'r--', 'b-', 'g-.', 'm:'};
t_eval = linspace(0, 14, 1401)';

fig2 = figure('Color', 'w', 'Position', [100, 100, 1000, 550]);
fprintf('\n--- Periodicity Residual ||u(t+omega) - u(t)||_E ---\n');

for m = 1:length(tols)
    opts_m = odeset('RelTol', tols(m), 'AbsTol', 1e-2 * tols(m));
    [tm, um] = ode45(@(t, u) eco_climatic_canonical_system(t, u), ...
        t_span, IC1, opts_m);
    u_t   = interp1(tm, um, t_eval, 'spline');
    u_tw  = interp1(tm, um, t_eval + omega, 'spline');
    err_m = vecnorm(u_tw - u_t, 2, 2);

    semilogy(t_eval, err_m, styles{m}, 'LineWidth', 3);
    hold on;
    fprintf('RelTol = %.0e | E_per(5) = %.3e | E_per(10) = %.3e | E_per(14) = %.3e\n', ...
        tols(m), err_m(t_eval == 5), err_m(t_eval == 10), err_m(end));
end

xlim([0 14]);
ylim([1e-12 1e2]);
legend('$\mathrm{RelTol}=10^{-4},\ \mathrm{AbsTol}=10^{-6}$', ...
    '$\mathrm{RelTol}=10^{-6},\ \mathrm{AbsTol}=10^{-8}$', ...
    '$\mathrm{RelTol}=10^{-8},\ \mathrm{AbsTol}=10^{-10}$', ...
    '$\mathrm{RelTol}=10^{-10},\ \mathrm{AbsTol}=10^{-12}$', ...
    'FontSize', 14, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Periodicity residual $\|\mathbf{u}(t+\omega)-\mathbf{u}(t)\|_E$', ...
    'FontSize', 14, 'Interpreter', 'latex');
title('Asymptotic $\omega$-Periodicity Convergence across ODE Tolerances', ...
    'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

file_path2_png = fullfile(folder_name, 'periodicity_error_diagnostic.png');
file_path2_eps = fullfile(folder_name, 'periodicity_error_diagnostic.eps');

set(fig2, 'PaperPositionMode', 'auto');
print(fig2, file_path2_png, '-dpng', '-r300'); 
print(fig2, file_path2_eps, '-depsc', '-r300'); 

fprintf('Done! Scenario A figures printed to "%s".\n', folder_name);

% =========================
% U system for Scenario A
% =========================
function dudt = eco_climatic_canonical_system(t, u)
    K = 10.0; 
    L = 8.0; 
    beta1 = 0.05; 
    beta2 = 0.2; 
    gamma1 = 0.4; 
    gamma2 = 0.3; 
    Phi1 = 0.5; 
    Phi2 = 0.4; 
    delta2 = 0.15; 
    rho2 = 0.25; 
    Q0 = 2.0; 
    tauC = 5.0; 
    lambda_val = 0.1; 
    alpha1 = 0.05; 
    alpha2 = 0.8;
    eta = 0.3; 
    h0 = 0.4; 
    gamma3 = 0.5; 
    KT = 5.0; n = 2;
    
    lambda_param = 1.0; 

    u1 = u(1); 
    u2 = u(2); 
    u3 = u(3); 
    u4 = u(4);
    
    % Parameter functions with omega = 1
    r1_t = 1.2 + 0.4 * sin(2 * pi * t);
    r2_t = 0.8 + 0.3 * cos(2 * pi * t);
    r3_t = 0.6 + 0.2 * sin(2 * pi * t);
    Tb_t = 15.0 + 5.0 * sin(2 * pi * t); 
    
    % System
    du1dt = -r1_t * u1 ...
        - lambda_param * (r3_t / K) * ...
        ((gamma1 - beta1) * K * u1^2) / (u2 * (u1 + Phi1)) ...
        + r1_t / K;

    du2dt = -(r2_t - delta2) * u2 ...
        + lambda_param * (r3_t / L) * ...
        ((rho2 * L * u2^2) / (u2 + beta2) ...
        - (gamma2 * L * u1 * u2) / (u1^2 + Phi2)) ...
        + r2_t / L;

    du3dt = -(1 / tauC) * u3 ...
        + lambda_param * ((eta * u4 * u3^n) / (u3^n + KT^n) ...
        - (lambda_val * u3) / (u1 * (1 + alpha1 * u3))) ...
        + lambda_param * Q0;

    du4dt = -alpha2 * u4 ...
        + lambda_param * ((h0 * u3^n) / (u3^n + KT^n) ...
        + gamma3 / (1 + u1 * K)) ...
        + lambda_param * alpha2 * Tb_t;

    dudt = [du1dt; du2dt; du3dt; du4dt];
end
