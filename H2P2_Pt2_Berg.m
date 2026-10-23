function flywheel_optimization()
    clc; clear; close all;

    % Conversion factor
    in2m = 0.0254;

    % Given constants (SI units: m, kg, s)
    params.rho = 7150;                  % Density (kg/m^3)
    params.r_s = 0.012;                 % Shaft radius (12 mm)
    params.r_h = 0.034;                 % Hub outer radius (34 mm)
    params.J_req = 0.1007659;           % Required moment of inertia (kg*m^2)
    params.min_tr = 1.5 * in2m;         % Min rim radial thickness (38.1 mm)
    params.min_t  = 1.0 * in2m;         % Min axial thickness (25.4 mm)
    params.max_ro = 8.0 * in2m;         % Max outer radius (203.2 mm)
    params.max_t  = 4.0 * in2m;         % Max axial thickness (101.6 mm)
    params.min_Ws = 0.010;              % Min feature size (10 mm)

    %% 1. Numerical Optimization using fmincon
    % Optimization variables x = [r_o, W_s, t, t_r]
    x0 = [0.150, 0.015, 0.030, 0.040];  % Initial guess
    
    % Bounds: [r_o, W_s, t, t_r]
    lb = [params.r_h + params.min_tr, params.min_Ws, params.min_t, params.min_tr];
    ub = [params.max_ro, 0.080, params.max_t, 0.150];

    % Call fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [x_opt, V_opt] = fmincon(@(x) obj_volume(x, params), x0, [], [], [], [], lb, ub, ...
                            @(x) nonlcon(x, params), options);

    r_o_opt = x_opt(1);
    W_s_opt = x_opt(2);
    t_opt   = x_opt(3);
    t_r_opt = x_opt(4);
    r_i_opt = r_o_opt - t_r_opt;
    mass_opt = V_opt * params.rho;

    % Display results
    fprintf('\n================ OPTIMAL FLYWHEEL DESIGN ================\n');
    fprintf('Outer Radius (r_o):         %.2f mm (%.3f in)\n', r_o_opt*1000, r_o_opt/in2m);
    fprintf('Spoke Width (W_s):          %.2f mm\n', W_s_opt*1000);
    fprintf('Axial Thickness (t):        %.2f mm (%.3f in)\n', t_opt*1000, t_opt/in2m);
    fprintf('Rim Radial Thickness (t_r): %.2f mm (%.3f in)\n', t_r_opt*1000, t_r_opt/in2m);
    fprintf('Rim Inner Radius (r_i):     %.2f mm\n', r_i_opt*1000);
    fprintf('Total Volume:               %.2f cm^3 (%.6f m^3)\n', V_opt*1e6, V_opt);
    fprintf('Total Mass:                 %.2f kg\n', mass_opt);
    fprintf('=========================================================\n\n');

    %% 2. Generate Contour Plot over (r_o, W_s) Space
    % Grid setup for r_o and W_s
    ro_grid = linspace(0.100, params.max_ro, 150); % m
    Ws_grid = linspace(0.008, 0.040, 150);          % m
    [RO, WS] = meshgrid(ro_grid, Ws_grid);

    % Calculate required thickness t and Volume for each grid point
    % assuming t_r = min_tr (38.1 mm) to minimize volume
    V_grid = zeros(size(RO));
    t_req_grid = zeros(size(RO));

    for i = 1:size(RO, 1)
        for j = 1:size(RO, 2)
            ro_val = RO(i,j);
            Ws_val = WS(i,j);
            ri_val = ro_val - params.min_tr;

            % Calculate J per unit thickness (J_prime) and V per unit thickness (V_prime)
            [V_prime, J_prime] = compute_unit_properties(ro_val, ri_val, Ws_val, params);
            
            % Required thickness to achieve J_req exactly
            t_req = params.J_req / J_prime;
            t_req_grid(i,j) = t_req;
            V_grid(i,j) = V_prime * t_req;
        end
    end

    % Plotting
    figure('Color', 'w', 'Position', [100, 100, 850, 650]);
    [C, h] = contourf(RO*1000, WS*1000, V_grid*1e6, 25, 'LineColor', 'none');
    colormap(jet);
    colorbar_h = colorbar;
    ylabel(colorbar_h, 'Total Volume (cm^3)', 'FontSize', 11, 'FontWeight', 'bold');
    hold on;

    % Contour lines for Volume
    [C2, h2] = contour(RO*1000, WS*1000, V_grid*1e6, 15, 'k--', 'LineWidth', 0.5);
    clabel(C2, h2, 'FontSize', 8, 'Color', 'k');

    % Constraint Boundaries
    % 1. Thickness constraint boundary: t_req = 1.0 inch (25.4 mm)
    contour(RO*1000, WS*1000, t_req_grid*1000, [25.4 25.4], 'r-', 'LineWidth', 2.5);

    % 2. Minimum feature size boundary (W_s = 10 mm)
    line([min(ro_grid)*1000, max(ro_grid)*1000], [10 10], 'Color', 'm', 'LineWidth', 2.5, 'LineStyle', '-.');

    % Mark Optimal Design Point
    plot(r_o_opt*1000, W_s_opt*1000, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 14, 'LineWidth', 1.5);

    % Annotations and Labels
    title('Flywheel Volume Optimization Contour Plot', 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Outer Radius r_o (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Spoke Width W_s (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    axis tight;

    legend({'Volume Contours', '', 'Active Boundary: t = 25.4 mm (1 in)', ...
            'Active Boundary: W_s = 10 mm', 'Optimal Design Point'}, ...
            'Location', 'northeast', 'FontSize', 9);
end

%% Helper Functions
function V = obj_volume(x, params)
    r_o = x(1); W_s = x(2); t = x(3); t_r = x(4);
    r_i = r_o - t_r;
    V_hub = pi * (params.r_h^2 - params.r_s^2) * t;
    V_rim = pi * (r_o^2 - r_i^2) * t;
    V_spokes = 4 * W_s * (r_i - params.r_h) * t;
    V = V_hub + V_rim + V_spokes;
end

function [c, ceq] = nonlcon(x, params)
    r_o = x(1); W_s = x(2); t = x(3); t_r = x(4);
    r_i = r_o - t_r;

    % Inertia calculation
    J_hub = 0.5 * params.rho * pi * t * (params.r_h^4 - params.r_s^4);
    J_rim = 0.5 * params.rho * pi * t * (r_o^4 - r_i^4);
    m_spoke = params.rho * W_s * t * (r_i - params.r_h);
    J_spoke = m_spoke * ((1/3)*(r_i^2 + params.r_h*r_i + params.r_h^2) + (1/12)*W_s^2);
    J_total = J_hub + J_rim + 4 * J_spoke;

    % Nonlinear constraints (c <= 0)
    c(1) = params.J_req - J_total;        % J >= J_req
    c(2) = params.r_h - r_i;              % r_i >= r_h
    ceq = [];
end

function [V_prime, J_prime] = compute_unit_properties(r_o, r_i, W_s, params)
    V_prime = pi * (params.r_h^2 - params.r_s^2) + pi * (r_o^2 - r_i^2) + 4 * W_s * (r_i - params.r_h);
    J_hub_prime = 0.5 * params.rho * pi * (params.r_h^4 - params.r_s^4);
    J_rim_prime = 0.5 * params.rho * pi * (r_o^4 - r_i^4);
    m_spoke_prime = params.rho * W_s * (r_i - params.r_h);
    J_spoke_prime = m_spoke_prime * ((1/3)*(r_i^2 + params.r_h*r_i + params.r_h^2) + (1/12)*W_s^2);
    J_prime = J_hub_prime + J_rim_prime + 4 * J_spoke_prime;
end