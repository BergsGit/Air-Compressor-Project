clc
clear
close all

%% Setup Parameters
target_J = 0.1008;       % Required J [kg*m^2]
b_fixed  = 0.0254;       % Fixed axial thickness [m] (1 inch)

params.rho   = 7150;     % Material density [kg/m^3]
params.rs    = 0.012;    % Shaft radius [m] (12 mm)
params.J_req = target_J; % Required moment of inertia [kg*m^2]
params.b_fixed = b_fixed;

%% Analytical Baseline Calculation (b = 1 in)
ro_solid = (((2 * target_J)/(params.rho * pi * b_fixed)) + (params.rs^4))^(1/4); 
mass_solid = 0.5 * params.rho * pi * b_fixed * (ro_solid^2 - params.rs^2);
fprintf('Analytical Check (b = %.2f mm / 1 in): r_o = %.2f mm, Mass = %.2f kg\n', ...
        b_fixed*1000, ro_solid*1000, mass_solid);

%% Run Optimization & Plotting
flywheel_optimization(params);

%% Optimization Function Definition
function flywheel_optimization(params)
    % Global variable to record fmincon history
    global x_history;
    
    b_target = params.b_fixed; % 0.0254 m (1 in)
    
    % Design Vector x = [r_o, b]
    min_ro = 0.020;  % Min outer radius [m]
    max_ro = 0.200;  % Max outer radius [m]
    
    % Constrain b to strictly equal 1 inch (0.0254 m)
    min_b = b_target;  
    max_b = b_target;  
    
    x0 = [0.100, b_target]; % Initial guess [r_o, b]
    x_history = x0;
    
    lowerbound = [min_ro, min_b];
    upperbound = [max_ro, max_b];
    
    % Optimization Options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'OutputFcn', @outfun);
                       
    % Run fmincon
    [x_opt, V_opt, exitflag, output] = fmincon(@(x) obj_volume(x, params), x0, ...
                             [], [], [], [], lowerbound, upperbound, ...
                             @(x) nonlcon(x, params), options);
                         
    ro_opt   = x_opt(1);
    b_opt    = x_opt(2);
    mass_opt = V_opt * params.rho;
    
    % Print Results
    fprintf('\n-----------------------------------\n');
    fprintf('OPTIMAL SOLID FLYWHEEL DESIGN RESULTS (b = 1")\n');
    fprintf('Iterations:           %d\n', output.iterations);
    fprintf('Function Evaluations: %d\n', output.funcCount);
    fprintf('-----------------------------------\n');
    fprintf('Optimal Outer Radius (r_o): %.2f mm\n', ro_opt*1000);
    fprintf('Fixed Thickness (b):        %.2f mm (%.2f in)\n', b_opt*1000, b_opt/0.0254);
    fprintf('Total Mass:                  %.2f kg\n', mass_opt);
    fprintf('-----------------------------------\n');

    %% Generate Contour Plot (Outer Radius vs Axial Thickness)
    ro_grid = linspace(min_ro, max_ro, 200);
    b_grid  = linspace(0.005, 0.050, 200); % Display window around 1 inch
    [RO, B] = meshgrid(ro_grid, b_grid);
    
    % Evaluate Volume and Inertia across the domain grid
    Vol_grid = pi * (RO.^2 - params.rs^2) .* B;
    J_grid   = 0.5 * params.rho * pi * B .* (RO.^4 - params.rs^4);
    
    % Plot Settings
    figure('Color', 'w', 'Position', [100, 100, 850, 650]);
    
    % 1. Volume Contours
    [C1, h1] = contour(RO*1000, B*1000, Vol_grid, 20, 'LineWidth', 1.2);
    cbar = colorbar;
    ylabel(cbar, 'Flywheel Volume (m^3)', 'FontSize', 11, 'FontWeight', 'bold');
    hold on;
    clabel(C1, h1, 'FontSize', 8, 'Color', 'k');
    
    % 2. Required Inertia Constraint Boundary (J = J_req)
    [C2, h2] = contour(RO*1000, B*1000, J_grid, [params.J_req, params.J_req], ...
                       'm', 'LineWidth', 1.5);
    
    % 3. Fixed Thickness Constraint Line (b = 1 inch)
    p_bline = yline(b_target*1000, 'm--', 'LineWidth', 1.5);
    
    % 4. Optimization Search Path
    p_steps   = plot(x_history(:,1)*1000, x_history(:,2)*1000, 'ro--', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    p_initial = plot(x0(1)*1000, x0(2)*1000, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
    
    % 5. Mark Optimum Point
    p_opt = plot(ro_opt*1000, b_opt*1000, 'o', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
    
    % Formatting
    xlabel('Outer Radius r_o (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Axial Thickness b (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    axis tight;
    
    legend([h1, h2, p_bline, p_steps, p_initial, p_opt], ...
           {'Volume Contours', 'Boundary (J_{req} = 0.1008 kg\cdotm^2)', ...
            'Boundary (Axial Thickness = 1 inch)', 'Path Iteration Steps', ...
            'Initial Guess', 'Optimal Design Point'}, ...
           'Location', 'Northeast');
end

%% Objective Function: Solid Flywheel Volume
function V = obj_volume(x, params)
    r_o = x(1);
    b   = x(2);
    V   = pi * (r_o^2 - params.rs^2) * b;
end

%% Constraint Function: Moment of Inertia Requirement
function [c, ceq] = nonlcon(x, params)
    r_o = x(1);
    b   = x(2);
    
    % Moment of inertia for solid flywheel with center hole
    J = 0.5 * params.rho * pi * b * (r_o^4 - params.rs^4);
    
    % Inequality constraint: J >= J_req  =>  J_req - J <= 0
    c = params.J_req - J;
    ceq = [];
end

%% Output Function to Track Progress
function stop = outfun(x, optimValues, state)
    stop = false;
    global x_history;
    if strcmp(state, 'iter')
        x_history = [x_history; x(1:2)];
    end
end