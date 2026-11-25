% =========================================================================
% MATLAB Script: Escaneo de sensibilidad y curvas para era-5487
% Description: Non-Autonomous Eco-Climatic Dynamics (Sensitivity Curves)
% Author: Anibal Coronel 
% =========================================================================

function sensibilidad_era5487()
    clear; clc; close all;

    % --- Parámetros base ---
    P = struct();
    P.K = 10.0; 
    P.L = 8.0; 
    P.beta1 = 0.05;
    P.beta2 = 0.2; 
    P.gamma1 = 0.4;
    P.gamma2 = 0.3; 
    P.Phi1 = 0.5; 
    P.Phi2 = 0.4; 
    P.delta2 = 0.15;
    P.rho2 = 0.25; 
    P.Q0 = 2.0; 
    P.tauC = 5.0;
    P.lam0 = 0.1; 
    P.alpha1 = 0.05;
    P.alpha2 = 0.8; 
    P.eta = 0.3; 
    P.h0 = 0.4; 
    P.gamma3 = 0.5; 
    P.KT = 5.0; 
    P.n = 2;
    P.canopy_sign = -1.0;  
    P.thermal = true;      
    P.Topt = 15.0; 
    P.Wth = 6.0; 

    % --- High  resolution ---
    N_res = 30;
    Q0_grid = linspace(0.5, 12.0, N_res);
    alpha2_grid = linspace(0.4, 2.0, N_res);
    gamma3_grid = linspace(0.0, 3.0, N_res);
    
    % Result  vectors
    res_Q0_Tmax = zeros(1, N_res);   
    res_Q0_ymean = zeros(1, N_res);
    res_a2_Tmean = zeros(1, N_res);  
    res_a2_ymean = zeros(1, N_res);
    res_g3_Tmax = zeros(1, N_res);   
    res_g3_ymin = zeros(1, N_res);
    
    fprintf('Runing ...\n');
    for i = 1:N_res
        % Scan Q0
        q = P; 
        q.Q0 = Q0_grid(i); 
        s = summary_orbit(q);
        res_Q0_Tmax(i) = s.Tmax; 
        res_Q0_ymean(i) = s.y;
        
        % Scan alpha2
        q = P; 
        q.alpha2 = alpha2_grid(i); 
        s = summary_orbit(q);
        res_a2_Tmean(i) = s.T; 
        res_a2_ymean(i) = s.y;
        
        % Scan gamma3
        q = P; 
        q.gamma3 = gamma3_grid(i); 
        s = summary_orbit(q);
        res_g3_Tmax(i) = s.Tmax; 
        res_g3_ymin(i) = s.ymin;
    end

    % --- MULTI-PANEL ---
    fig = figure('Color', 'w', 'Position', [100, 100, 1200, 450]);
    
    % Panel 1: Emisiones vs Polinizadores e Impacto Térmico
    subplot(1, 3, 1);
    yyaxis left
    h1 = plot(Q0_grid, res_Q0_ymean, 'r-', 'LineWidth', 3);
    ylabel('Pollinator Mean Density $\bar{y}$', ...
        'Interpreter', 'latex');
    set(gca, 'YColor', 'r'); 
    ylim([6.595, 6.615]);

    yyaxis right
    h2 = plot(Q0_grid, res_Q0_Tmax, 'k--', 'LineWidth', 3);
    ylabel('Max Temperature $T_{\max}$ ($^\circ$C)', ...
        'Interpreter', 'latex');
    set(gca, 'YColor', 'k'); 
    ylim([15.1, 15.9]);
    
    xlim([0.5, 12]);
    xlabel('GHG Emission Rate $Q_0$', 'Interpreter', 'latex');
    title('Sensitivity on Anthropogenic Emissions', ...
        'Interpreter', 'latex');
    legend([h1, h2], {'Mean Density $\bar{y}$', ...
        'Max Temperature $T_{\max}$'}, ...
           'FontSize', 10, 'Interpreter', 'latex', ...
           'Location', 'southeast');
    grid on; set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 14, 'LineWidth', 1.5);
    
    % Panel 2: Disipación Térmica
    subplot(1, 3, 2);
    plot(alpha2_grid, res_a2_Tmean, 'b-', 'LineWidth', 3);
    xlabel('Thermal Conductivity $\alpha_2$', ...
        'Interpreter', 'latex');
    ylabel('Mean Temperature Elevation $\bar{T}$ ($^\circ$C)', ...
        'Interpreter', 'latex');
    title('Sensitivity on Thermal Conductivity', ...
        'Interpreter', 'latex');
    legend('Mean Temperature $\bar{T}$', 'FontSize', ...
        14, 'Interpreter', 'latex', 'Location', 'northwest');
    xlim([0.4, 2.0]);        
    ylim([14.4, 14.95]);     
    set(gca, 'XTick', 0.4:0.4:2.0);
    grid on; set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 14, 'LineWidth', 1.5);
    set(gca, 'YColor', 'k');
    
    % Panel 3: Mitigación por Sombra Forestal (Canopy)
    subplot(1, 3, 3);
    yyaxis left
    h3 = plot(gamma3_grid, res_g3_ymin, 'g-', 'LineWidth', 3);
    ylabel('Pollinator Min Density $y_{\min}$', ...
        'Interpreter', 'latex');
    set(gca, 'YColor', [0 0.5 0]); 
    ylim([6.05, 6.6]);
    
    yyaxis right
    h4 = plot(gamma3_grid, res_g3_Tmax, 'k--', 'LineWidth', 3);
    ylabel('Max Temperature $T_{\max}$ ($^\circ$C)', ...
        'Interpreter', 'latex');
    set(gca, 'YColor', 'k'); 
    ylim([13.3, 15.7]);
    xlim([0.0, 3.0]);       
    set(gca, 'XTick', 0.0:0.5:3.0);
    
    xlabel('Microclimatic Buffer Rate $\gamma_3$', ...
        'Interpreter', 'latex');
    title('Sensitivity on Canopy Protection', ...
        'Interpreter', 'latex');
    legend([h3, h4], {'Min Density $y_{\min}$', ...
        'Max Temperature $T_{\max}$'}, ...
           'FontSize', 14, 'Interpreter', 'latex', ...
           'Location', 'southwest');
    grid on; set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 14, 'LineWidth', 1.5);

    % --- Exportación a carpeta FIGURASERA ---
    folder_name = 'FIGURASERA';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    print(fig, fullfile(folder_name, 'sensitivity_curves_scan.png'), ...
        '-dpng', '-r300');
    print(fig, fullfile(folder_name, 'sensitivity_curves_scan.eps'), ...
        '-depsc', '-r300');
    fprintf('Done! figures printed to "%s".\n', folder_name);
end

% =========================================================================
% SUBFUNCIONES INTERNAS DE DINÁMICA ECO-CLIMÁTICA
% =========================================================================

function out = r1(t)
    out = 1.2 + 0.4 * sin(2 * pi * t);
end

function out = r2(t)
    out = 0.8 + 0.3 * cos(2 * pi * t);
end

function out = r3(t)
    out = 0.6 + 0.2 * sin(2 * pi * t);
end

function out = Tb(t)
    out = 15.0 + 5.0 * sin(2 * pi * t);
end

function out = psi_func(T, p)
    if ~p.thermal
        out = 1.0;
    else
        out = exp(-((T - p.Topt) / p.Wth)^2);
    end
end

function dz = rhs(t, z, p)
    x = z(1); y = z(2); C = z(3); T = z(4);
    hill = (C^p.n) / (C^p.n + p.KT^p.n);
    s = psi_func(T, p);
    
    dx = r1(t)*x*(1 - x/p.K) + (p.gamma1 - ...
        p.beta1)*r3(t)*s*x*y / (1 + p.Phi1*x);
    dy = (r2(t)*s*y*(1 - y/p.L) + ...
        p.gamma2*r3(t)*s*(x^2)*y / (1 + p.Phi2*(x^2)) ...
          - p.delta2*y - p.rho2*r3(t)*(y^2) / (1 + p.beta2*y));
    dC = p.Q0 - C/p.tauC - ...
        p.lam0*x*C / (1 + p.alpha1*C) + p.eta*T*hill;
    dT = p.h0*hill - ...
        p.alpha2*(T - Tb(t)) + p.canopy_sign*p.gamma3*x / (x + p.K);
    
    dz = [dx; dy; dC; dT];
end

function z = rk4_integrator(z0, t0, t1, p, h)
    if nargin < 5, h = 1e-3; end
    z = double(z0); t = t0;
    N = round((t1 - t0) / h);
    for k = 1:N
        k1 = rhs(t, z, p);
        k2 = rhs(t + h/2, z + (h/2)*k1, p);
        k3 = rhs(t + h/2, z + (h/2)*k2, p);
        k4 = rhs(t + h, z + h*k3, p);
        z = z + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
    end
end

function traj = orbit_trajectory(z0, t0, p, h)
    if nargin < 4, h = 1e-3; end
    z = double(z0); t = t0;
    N = round(1.0 / h); 
    traj = zeros(N + 1, 4);
    traj(1, :) = z';
    for k = 1:N
        k1 = rhs(t, z, p);
        k2 = rhs(t + h/2, z + (h/2)*k1, p);
        k3 = rhs(t + h/2, z + (h/2)*k2, p);
        k4 = rhs(t + h, z + h*k3, p);
        z = z + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
        traj(k + 1, :) = z';
    end
end

function [zT, traj, mult] = periodic_orbit(p)
    z0 = [5.0; 4.0; 5.0; 15.0];
    transient = 60.0; h = 2e-3;
    zT = rk4_integrator(z0, 0.0, transient, p, h);
    
    for iter = 1:4
        F = rk4_integrator(zT, 0.0, 1.0, p, h) - zT;
        J = zeros(4, 4);
        for j = 1:4
            e = zeros(4, 1);
            e(j) = 1e-6 * max(1.0, abs(zT(j)));
            J(:, j) = (rk4_integrator(zT + e, 0.0, 1.0, p, h) - ...
                rk4_integrator(zT - e, 0.0, 1.0, p, h)) / (2 * e(j));
        end
        M = J; 
        dz = (M - eye(4)) \ (-F);
        zT = zT + dz;
        if norm(dz) < 1e-10
            break;
        end
    end
    mult = eig(M);
    traj = orbit_trajectory(zT, 0.0, p, h);
end

function s = summary_orbit(p)
    [~, traj, mult] = periodic_orbit(p);
    m = mean(traj, 1);
    s = struct();
    s.x = m(1); s.y = m(2); s.C = m(3); s.T = m(4);
    s.Tmax = max(traj(:, 4));
    s.ymin = min(traj(:, 2));
    s.rho = max(abs(mult));
end
