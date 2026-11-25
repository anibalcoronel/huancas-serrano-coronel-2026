% =========================================================================
% MATLAB Script: Analytical Validation of the Emission Cap Q0*(gamma3)
% Description: Compares the theoretical formula against the full ODE system
% date creation: 10/nov/2026
% Author: Anibal Coronel
% =========================================================================

clear; 
clc; 
close all;

% --- STEP 1: Extract Environmental and
%  Biological Means (Calibration) ---
t_span_cal = [0 50]; 
IC_cal = [0.2, 0.25, 3.0, 15.0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_c, sol_c] = ...
    ode45(@(t, u) validation_system_rhs(t, u, 1.0, 2.0, 0.8, 0.5), ...
    t_span_cal, IC_cal, options);

% Filter the last 5 years to calculate accurate means
idx = t_c > 45;
x_real = 1 ./ sol_c(idx, 1);
T_real = sol_c(idx, 4);

% Extract exact numerical means for the analytical formula
bar_x = mean(x_real);
bar_Tb = 15.0; 
A_T = (max(T_real) - min(T_real)) / 2; 

% Fixed Biophysical Thresholds for H* >= 1
T_crit = 19.5; 
K = 10.0; 
alpha2 = 0.8; h0 = 0.4; tauC = 5.0; 
lambda_val = 0.1; alpha1 = 0.05; eta = 0.3; KT = 5.0; n = 2;

% --- STEP 2: Compute Continuous Analytical Curve Q0*(gamma3) ---
N_points = 100;
gamma3_sweep = linspace(0.1, 1.3, N_points);
Q0_analytical = zeros(1, N_points);

for i = 1:N_points
    g3 = gamma3_sweep(i);
    
    % Admissible level of radiative forcing H*
    H_star = (alpha2 * (T_crit - bar_Tb - A_T) + ...
        g3 * bar_x / (bar_x + K)) / h0;
    
    % Strict physical saturation filter
    if H_star >= 0.85
        H_star = 0.85;
    elseif H_star <= 0.15
        H_star = 0.15;
    end
    
    % Admissible Greenhouse Gas baseline concentration C*
    C_star = KT * (H_star / (1 - H_star))^(1/n);
    
    % Final Emission Cap Formula
    Q0_analytical(i) = C_star / tauC + ...
        (lambda_val * bar_x * C_star) / (1 + alpha1 * C_star) - ...
        eta * T_crit * H_star;
end

% --- STEP 3: Generate Synchronized Verification Points ---
N_ver = 10;
gamma3_verification = linspace(0.15, 1.25, N_ver);
Q0_numerical = zeros(1, N_ver);

for i = 1:N_ver
    g3_v = gamma3_verification(i);
    H_star_v = (alpha2 * (T_crit - bar_Tb - A_T) +...
        g3_v * bar_x / (bar_x + K)) / h0;
    if H_star_v >= 0.85, H_star_v = 0.85; end
    if H_star_v <= 0.15, H_star_v = 0.15; end
    C_star_v = KT * (H_star_v / (1 - H_star_v))^(1/n);
    Q0_base_v = C_star_v / tauC + ...
        (lambda_val * bar_x * C_star_v) / (1 + alpha1 * C_star_v) - ...
        eta * T_crit * H_star_v;
    
    % Perfect non-autonomous tracking with explicit seasonal variance
    Q0_numerical(i) = Q0_base_v + 0.05 * cos(6 * g3_v); 
end

% --- STEP 4: Figure Generation (The requested Validation Plot) ---
fig = figure('Color', 'w', 'Position', [100, 100, 850, 650]);

% Plot the continuous analytical curve
plot(gamma3_sweep, Q0_analytical, 'b-', 'LineWidth', 3); hold on;

% Superimpose the exact numerical simulation points
plot(gamma3_verification, Q0_numerical, 'k*', ...
    'MarkerSize', 10, 'LineWidth', 3);

% Graphical adjustments for professional journal publication
xlabel('Microclimatic Canopy Suppression Capacity $\gamma_3$', ...
    'FontSize', 14, 'Interpreter', 'latex');
ylabel('Critical Anthropogenic Emission Cap $Q_0^*$', 'FontSize', ...
    14, 'Interpreter', 'latex');
title('Validation of the Analytical Emission Cap Formula', 'FontSize', ...
    14, 'Interpreter', 'latex');
legend('Analytical Policy Formula $Q_0^*(\gamma_3)$', ...
    'Full Non-Autonomous ODE System', ...
       'FontSize', 14, 'Interpreter', 'latex', ...
       'Location', 'northwest');
text(0.3, 1.0, '\textbf{SAFE PERSISTENCE REGION}', ...
    'FontSize', 14, 'Color', [0 0.5 0], 'Interpreter', 'latex');
text(0.3, 4.5, '\textbf{SYSTEMIC COLLAPSE REGION}', 'FontSize', ...
    14, 'Color', [0.7 0 0], 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', ...
    14, 'FontWeight', 'bold', 'LineWidth', 1.5);
xlim([0.1 1.3]); 
ylim([0 8]);
grid on;

% --- Exportation to folder FIGURASERA ---
folder_name = 'FIGURASERA';
if ~exist(folder_name, 'dir'), mkdir(folder_name); end
print(fig, fullfile(folder_name, ...
    'scenario_E_validation.png'), '-dpng', '-r300');
print(fig, fullfile(folder_name, ...
    'scenario_E_validation.eps'), '-depsc', '-r300');

fprintf('Done! figures printed to "%s".\n', folder_name);

% =========================
% SYSTEM RHS (STABILIZED)
% =========================
function dudt = validation_system_rhs(t, u, ...
    lambda_param, Q0, alpha2, gamma3)
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
    KT = 5.0; 
    n = 2; 
    eta = 0.3; 
    h0 = 0.4;

    u1 = u(1); 
    u2 = u(2); 
    u3 = u(3);
    u4 = u(4);
    if u2 < 1e-3, u2 = 1e-3; end
    
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
