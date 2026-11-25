% =========================================================================
% MATLAB Script: ECO-Climatic plant--pollinator (Canonical Scenario C) 
% date creation: 5/oct/2026
% date modification: 8/oct/2026
% Author: Anibal Coronel
% =========================================================================

clear; 
clc; 
close all;

% --- Time interval for simulation ---
t_span = [0 25]; 

% Identical Initial Conditions for both cases in Scenario C (1/x0, 1/y0, C0, T0)
IC = [1/5.0, 1/4.0, 3.0, 15.0]; 

% --- Numerical solution with ODE45 ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t1, sol1] = ode45(@(t, u) eco_climatic_canonical_uno(t, u), ...
    t_span, IC, options);
[t2, sol2] = ode45(@(t, u) eco_climatic_canonical_dos(t, u), ...
    t_span, IC, options);

% --- Figures generation ---
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 750]);

subplot(2, 2, 1);
plot(t1, sol1(:, 1), 'r-', t2, sol2(:, 1), 'b:', 'LineWidth', 3);
legend('$\lambda = 0.5 < \lambda_c$ (Sustainable)', ...
    '$\lambda = 2.5 > \lambda_c$ (Over-forced)', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'east');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_1(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Plants $u_1(t) = 1/x(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 2);
plot(t1, sol1(:, 2), 'r-', t2, sol2(:, 2), 'b:', 'LineWidth', 3);
legend('$\lambda = 0.5 < \lambda_c$ (Sustainable)', ...
    '$\lambda = 2.5 > \lambda_c$ (Over-forced)', 'FontSize', 12, ...
    'Interpreter', 'latex', 'Location', 'east');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_2(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Transformed Pollinators $u_2(t) = 1/y(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 3);
plot(t1, sol1(:, 3), 'r-', t2, sol2(:, 3), 'b:', 'LineWidth', 3);
legend('$\lambda = 0.5 < \lambda_c$ (Sustainable)', ...
    '$\lambda = 2.5 > \lambda_c$ (Over-forced)', ...
    'FontSize', 12, 'Interpreter', 'latex', 'Location', 'east');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_3(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Greenhouse Gases $u_3(t) = C(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

subplot(2, 2, 4);
plot(t1, sol1(:, 4), 'r-', t2, sol2(:, 4), 'b:', 'LineWidth', 3);
legend('$\lambda = 0.5 < \lambda_c$ (Sustainable)', ...
    '$\lambda = 2.5 > \lambda_c$ (Over-forced)', 'FontSize', ...
    12, 'Interpreter', 'latex', 'Location', 'east');
xlabel('Time $t$ (years)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$u_4(t)$', 'FontSize', 14, 'Interpreter', 'latex');
title('Temperature Elevation $u_4(t) = T(t)$', 'FontSize', ...
    14, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on; 

% --- Exportation to folder ---
folder_name = 'FIGURASERA';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

file_path_png = fullfile(folder_name, 'scenario_C_bifurcation.png');
file_path_eps = fullfile(folder_name, 'scenario_C_bifurcation.eps');

set(fig, 'PaperPositionMode', 'auto');
print(fig, file_path_png, '-dpng', '-r300'); 
print(fig, file_path_eps, '-depsc', '-r300'); 

fprintf('Done! Scenario C figures printed to "%s".\n', folder_name);

% ================================================
% U SYSTEM ONE: SUSTAINABLE REGIME (lambda = 0.5)
% ================================================
function dudt = eco_climatic_canonical_uno(t, u)
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
    KT = 5.0; 
    n = 2;
    
    % Sustainable framework mapping under the threshold
    lambda_param = 0.5; 

    u1 = u(1); 
    u2 = u(2); 
    u3 = u(3); 
    u4 = u(4);
    
    r1_t = 1.2 + 0.4 * sin(2 * pi * t);
    r2_t = 0.8 + 0.3 * cos(2 * pi * t);
    r3_t = 0.6 + 0.2 * sin(2 * pi * t);
    Tb_t = 15.0 + 5.0 * sin(2 * pi * t);
        
    du1dt = -r1_t * u1 ...
        - lambda_param * (r3_t / K) * ...
        ((gamma1 - beta1) * K * u1^2) / (u2 * (u1 + Phi1)) ...
        + r1_t / K;

    du2dt = -(r2_t - delta2) * u2 ...
        + lambda_param * (r3_t / L) * ((rho2 * L * u2^2) / (u2 + beta2) ...
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

% ================================================
% U SYSTEM TWO: SUSTAINABLE REGIME (lambda = 2.5)
% ================================================

function dudt = eco_climatic_canonical_dos(t, u)
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
    KT = 5.0; 
    n = 2;
    
    % Boundary bifurcation forcing parameter
    lambda_param = 2.5; 

    u1 = u(1); 
    u2 = u(2); 
    u3 = u(3); 
    u4 = u(4);
    
    r1_t = 1.2 + 0.4 * sin(2 * pi * t);
    r2_t = 0.8 + 0.3 * cos(2 * pi * t);
    r3_t = 0.6 + 0.2 * sin(2 * pi * t);
    Tb_t = 15.0 + 5.0 * sin(2 * pi * t);
 
    du1dt = -r1_t * u1 ...
        - lambda_param * (r3_t / K) ...
        * ((gamma1 - beta1) * K * u1^2) / (u2 * (u1 + Phi1)) ...
        + r1_t / K;

    du2dt = -(r2_t - delta2) * u2 ...
        + lambda_param * (r3_t / L) ...
        * ((rho2 * L * u2^2) / (u2 + beta2) ...
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
