% =========================================================================
% MATLAB Script: ECO-Climatic plant--pollinator (Canonical Scenario B) 
% date creation: 5/oct/2026
% date modification: 8/oct/2026
% Author: Anibal Coronel
% =========================================================================

clear; 
clc; 
close all;

% --- Time interval for simulation ---
t_span = [0,1.5]; 

% Solution 1: Transformed Low initial conditions (1/x0, 1/y0, C0, T0)
% Real states: x0 = 0.8, y0 = 0.5, C0 = 0.5, T0 = 10.0
IC1 = [1/0.8, 1/0.5, 0.5, 10.0]; 

% Solution 2: Transformed High initial conditions (1/x0, 1/y0, C0, T0)
% Real states: x0 = 8.0, y0 = 6.0, C0 = 4.0, T0 = 18.0
IC2 = [1/8.0, 1/6.0, 4.0, 18.0]; 

% --- Numerical solution with ODE15s ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t1, sol1] = ode15s(@(t, u) eco_climatic_superlinear_system(t, u), ...
    t_span, IC1, options);
[t2, sol2] = ode15s(@(t, u) eco_climatic_superlinear_system(t, u), ...
    t_span, IC2, options);

% --- Figures generation ---
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

subplot(2, 2, 1);
plot(t1, sol1(:, 1), 'r-', t2, sol2(:, 1), 'b:', 'LineWidth', 3);
legend('Low-Amplitude Orbit', 'High-Amplitude Orbit', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_1(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Plants $u_1(t) = 1/x(t)$', 'FontSize', 14, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 2);
plot(t1, sol1(:, 2), 'r-', t2, sol2(:, 2), 'b:', 'LineWidth', 3);
legend('Low-Amplitude Orbit', 'High-Amplitude Orbit', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_2(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Pollinators $u_2(t) = 1/y(t)$', ...
    'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 3);
plot(t1, sol1(:, 3), 'r-', t2, sol2(:, 3), 'b:', 'LineWidth', 3);
legend('Low-Amplitude Orbit', 'High-Amplitude Orbit', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_3(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Greenhouse Gases $u_3(t) = C(t)$', 'FontSize', 14, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 4);
plot(t1, sol1(:, 4), 'r-', t2, sol2(:, 4), 'b:', 'LineWidth', 3);
legend('Low-Amplitude Orbit', 'High-Amplitude Orbit', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_4(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Temperature Elevation $u_4(t) = T(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

% --- Exportation to folder ---
folder_name = 'FIGURASERA';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

file_path_png = fullfile(folder_name, 'scenario_B_multiplicity.png');
file_path_eps = fullfile(folder_name, 'scenario_B_multiplicity.eps');

set(fig, 'PaperPositionMode', 'auto');
print(fig, file_path_png, '-dpng', '-r300'); 
print(fig, file_path_eps, '-depsc', '-r300'); 

fprintf('Done! Scenario B figures printed to "%s".\n', folder_name);

% =========================
% U system for Scenario B
% =========================
function dudt = eco_climatic_superlinear_system(t, u)
    % Base parameters
    K = 10.0; 
    L = 8.0; 
    beta1 = 0.05; 
    beta2 = 0.2; 
    gamma1 = 0.4; 
    gamma2 = 0.3; 
    delta2 = 0.15; 
    rho2 = 0.25; 
    Q0 = 2.0; 
    tauC = 5.0; 
    lambda_val = 0.1; 
    alpha2 = 0.8;
    eta = 0.3; 
    h0 = 0.4; 
    gamma3 = 0.5; 
    n = 2;
    
    % Small scaling value to balance 
    % topological multiplicity (Assertion b)
    lambda_param = 0.05; 

    u1 = u(1); 
    u2 = u(2); 
    u3 = u(3); 
    u4 = u(4);
    
    % Seasonal forcing functions
    r1_t = 1.2 + 0.4 * sin(2 * pi * t);
    r2_t = 0.8 + 0.3 * cos(2 * pi * t);
    r3_t = 0.6 + 0.2 * sin(2 * pi * t);
    Tb_t = 15.0 + 5.0 * sin(2 * pi * t);
    
    % System
    du1dt = -r1_t * u1 ...
        - lambda_param * (r3_t / K) * ...
        (-(gamma1 - beta1) * K * (u1 / u2)) ...
        + r1_t / K;

    du2dt = -(r2_t - delta2) * u2 ...
        + lambda_param * (r3_t / L) * (rho2 * L * u2^2 ...
        - gamma2 * L * u1 * u2) ...
        + r2_t / L;

    du3dt = -(1 / tauC) * u3 ...
        + lambda_param * (eta * u4 * u3^n ...
        - lambda_val * (u3 / u1)) ...
        + lambda_param * Q0;

    du4dt = -alpha2 * u4 ...
        + lambda_param * (h0 * u3^n ...
        + gamma3 * u1 * K) ...
        + lambda_param * alpha2 * Tb_t;

    
    dudt = [du1dt; du2dt; du3dt; du4dt];
end
