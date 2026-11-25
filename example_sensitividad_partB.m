% =========================================================================
% MATLAB Script: Quantitative Sensitivity Scan for Critical Threshold \lambda_c
% Description: Generates smooth multi-panel sensitivity curves for Q0, \alpha_2, \gamma_3
% Modified under Option A (h0 = 6.0) + Bisection Search (Step-free curves)
% date creation: 10/nov/2026
% Author: Anibal Coronel
% =========================================================================

clear; 
clc; 
close all;

% --- Parametric Scan Configuration ---
N_points = 50; % Resolution of the parameter sweep
lambda_bounds = [0.5, 4.0]; 
% Search interval [lambda_min, lambda_max] for bisection
t_span = [0 25]; 
IC = [0.2, 0.25, 3.0, 15.0]; 
% Standard Scenario C Initial Conditions
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

% Baseline parameters
Q0_base = 2.0; alpha2_base = 0.8; gamma3_base = 0.5;

% --- SCAN 1: Varying Greenhouse Gas Input (Q0) ---
Q0_sweep = linspace(1.0, 5.0, N_points);
lambda_c_Q0 = zeros(1, N_points);

for i = 1:N_points
    Q0_current = Q0_sweep(i);
    lambda_c_Q0(i) = ...
        find_bifurcation_threshold(t_span, IC, lambda_bounds, ...
        Q0_current, alpha2_base, gamma3_base, options);
end

% --- SCAN 2: Varying Thermodynamic Cooling (alpha2) ---
alpha2_sweep = linspace(0.4, 1.5, N_points);
lambda_c_alpha2 = zeros(1, N_points);

for i = 1:N_points
    alpha2_current = alpha2_sweep(i);
    lambda_c_alpha2(i) = ...
        find_bifurcation_threshold(t_span, IC, lambda_bounds, ...
        Q0_base, alpha2_current, gamma3_base, options);
end

% --- SCAN 3: Varying Microclimatic Canopy Buffering (gamma3) ---
gamma3_sweep = linspace(0.0, 1.5, N_points);
lambda_c_gamma3 = zeros(1, N_points);

for i = 1:N_points
    gamma3_current = gamma3_sweep(i);
    lambda_c_gamma3(i) = ...
        find_bifurcation_threshold(t_span, IC, lambda_bounds, ...
        Q0_base, alpha2_base, gamma3_current, options);
end

% --- FIGURES GENERATION (Multi-panel Figure 8b) ---
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(Q0_sweep, lambda_c_Q0, 'r-', 'LineWidth', 3);
xlabel('GHG Emission Rate $Q_0$', 'FontSize', 14, ...
    'Interpreter', 'latex');
ylabel('Critical Threshold $\lambda_c$', 'FontSize', 14, ...
    'Interpreter', 'latex');
title('Sensitivity Scan on $Q_0$', 'FontSize', 14, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

subplot(1, 3, 2);
plot(alpha2_sweep, lambda_c_alpha2, 'b-', 'LineWidth', 3);
xlabel('Thermal Conductivity $\alpha_2$', 'FontSize', 14, ...
    'Interpreter', 'latex');
ylabel('Critical Threshold $\lambda_c$', 'FontSize', 14, ...
    'Interpreter', 'latex');
title('Sensitivity Scan on $\alpha_2$', 'FontSize', 14, ...
    'Interpreter', 'latex');
xlim([0.4 1.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

subplot(1, 3, 3);
plot(gamma3_sweep, lambda_c_gamma3, 'g-', 'LineWidth', 3);
xlabel('Microclimatic Buffer Rate $\gamma_3$', 'FontSize', 14, ...
    'Interpreter', 'latex');
ylabel('Critical Threshold $\lambda_c$', 'FontSize', 14, ...
    'Interpreter', 'latex');
title('Sensitivity Scan on $\gamma_3$', 'FontSize', 14, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
    'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

% --- Exportation to folder ---
folder_name = 'FIGURASERA';
if ~exist(folder_name, 'dir'), mkdir(folder_name); end
print(fig, fullfile(folder_name, 'scenario_D_sensitivity.png'), ...
    '-dpng', '-r300');
print(fig, fullfile(folder_name, 'scenario_D_sensitivity.eps'), ...
    '-depsc', '-r300');

fprintf('Done! figures printed to "%s".\n', folder_name);

% =========================================================================
% AUXILIARY FUNCTION: BISECTION BIFURCATION SEARCH (SMOOTH CURVES)
% =========================================================================
function lambda_c = find_bifurcation_threshold(t_span, IC, ...
    lambda_bounds, Q0, alpha2, gamma3, options)

    lam_low = lambda_bounds(1);
    lam_high = lambda_bounds(2);
    
    % 16 iterations give precision ~ (4.0 - 0.5) / 2^16 = 5.3e-5
    for iter = 1:16
        lam_mid = 0.5 * (lam_low + lam_high);
        [t, sol] = ...
            ode45(@(t, u) scan_system_rhs(t, u, ...
            lam_mid, Q0, alpha2, gamma3), t_span, IC, options);
        
        % Evaluate collapse a(t > 2.0) with threshold u2 > 1.8
        if any(sol(t > 2.0, 2) > 1.8) || any(isnan(sol(:, 2)))
            lam_high = lam_mid; % Instability triggered -> lower upper bound
        else
            lam_low = lam_mid;  % Stable persistence -> raise lower bound
        end
    end
    lambda_c = 0.5 * (lam_low + lam_high);
end

% =========================================================================
% LOCAL SYSTEM EQUATIONS FOR SCANNING (OPTION A: h0 = 6.0)
% =========================================================================
function dudt = scan_system_rhs(t, u, lambda_param, Q0, alpha2, gamma3)
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
    tauC = 5.0; 
    lambda_val = 0.1; 
    alpha1 = 0.05; 
    KT = 5.0; n = 2;
    eta = 0.3; 
    h0 = 6.0;

    u1 = u(1); u2 = u(2); u3 = u(3); u4 = u(4);
    
    % Thermal curves (H8)
    T_opt1 = 15.0;
    T_opt2 = 18.0;
    sigma_T1 = 12.0; 
    sigma_T2 = 10.0; 
    psi1 = 0.05 + 0.95 * exp(-((u4 - T_opt1)^2) / (2 * sigma_T1^2));
    psi2 = 0.05 + 0.95 * exp(-((u4 - T_opt2)^2) / (2 * sigma_T2^2));
    
    r1_t = 1.2 + 0.4 * sin(2 * pi * t);
    r2_t = 0.8 + 0.3 * cos(2 * pi * t);
    r3_t = 0.6 + 0.2 * sin(2 * pi * t);
    Tb_t = 15.0 + 5.0 * sin(2 * pi * t); 
    
    du1dt = -r1_t * psi1 * u1 ...
        - lambda_param * (r3_t * psi2 / K) ...
        * ((gamma1 - beta1) * K * u1^2) / (u2 * (u1 + Phi1)) ...
        + r1_t * psi1 / K;

    du2dt = -(r2_t * psi2 - delta2) * u2 ...
        + lambda_param * (r3_t * psi2 / L) ...
        * ((rho2 * L * u2^2) / (u2 + beta2) ...
        - (gamma2 * L * u1 * u2) / (u1^2 + Phi2)) ...
        + r2_t * psi2 / L;

    du3dt = -(1 / tauC) * u3 ...
        + lambda_param * ((eta * u4 * u3^n) / (u3^n + KT^n) ...
        - (lambda_val * u3) / (u1 * (1 + alpha1 * u3))) ...
        + lambda_param * Q0;

    du4dt = -alpha2 * u4 ...
        + lambda_param * ((h0 * u3^n) / (u3^n + KT^n) ...
        - gamma3 / (1 + u1 * K)) ...
        + lambda_param * alpha2 * Tb_t;
    
    dudt = [du1dt; du2dt; du3dt; du4dt];
end
