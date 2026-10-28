%% Homework 2, Problem 2, Part 1
% Goal: to optimize the design of a solid flywheel to minimize mass
% while also meeting design requirements and parameters. 

% Clear workspace, command window, and any plots in MATLAB
clc
clear
close all

%% Setup Parameters
target_J = 0.1008; % Required J [kg*m^2]
b  = 0.0254; % Axial thickness [m] (1 inch)

% Define constants in a structure
params.rho   = 7150;           % Material density [kg/m^3]
params.rs    = 0.012;          % Shaft radius [m] (12 mm)
params.J_req = target_J;       % Required moment of inertia [kg*m^2]
params.b_fixed = b;            % Axial thickness [m]
 
%% Solve solid flywheel
% Required outer radius value
ro_solid = (((2 * target_J)/(params.rho * pi * b)) + (params.rs^4))^(1/4);
% Require mass value
mass_solid = 0.5 * params.rho * pi * b * (ro_solid^2 - params.rs^2);

% Call the optimization function
flywheel_optimization(params);

%% Optimization Function Definition
function flywheel_optimization(params)
    % Global array to store iteration history
    global x_steps;

    % Unpack structure
    b_target = params.b_fixed; % 0.0254 m (1 in)
    
    % Design Constraints
    min_ro = 0.020;  % Min outer radius
    max_ro = 0.2032;  % Max outer radius (8 inches)
    
    % Constrain b to equal 1 inch (0.0254 m)
    min_b = b_target;  
    max_b = b_target;  
    
    % Initial guess for fmincon
    x0 = [0.100, b_target]; % Initial guess [r_o, b]
    x_steps = x0; % Initialize itertation history with starting point
    
    % Defining upper and lower bounds
    lowerbound = [min_ro, min_b];
    upperbound = [max_ro, max_b];
    
    % Optimization Options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'OutputFcn', @outfun);
                       
    % Run fmincon
    [x_opt, V_opt, exitflag, output] = fmincon(@(x) obj_volume(x, params), x0, ...
                             [], [], [], [], lowerbound, upperbound, ...
                             @(x) nonlcon(x, params), options);
    % Store optimized results                   
    ro_opt   = x_opt(1); % Optimized outer radius [m]
    b_opt    = x_opt(2); % Optimized axial thickness [m]
    mass_opt = V_opt * params.rho; % Optimized mass [kg]
    
    % Print Results
    fprintf('Optimal Solid Flywheel Design Results (b = 1")\n');
    fprintf('Iterations: %d\n', output.iterations);
    fprintf('Function Evaluations: %d\n', output.funcCount);
    fprintf('Total Mass: %.2f kg\n', mass_opt);


    %% Generate Contour Plot over (r_o, b) Space
    ro_grid = linspace(min_ro, max_ro, 200);  %Outer Radius values
    b_grid = linspace(0.005, 0.050, 200); % Axial thickness values
    [RO, B] = meshgrid(ro_grid, b_grid);
    
    % Evaluate Volume and Inertia across the domain grid
    Vol_grid = pi * (RO.^2 - params.rs^2) .* B; % Volume grid
    J_grid = 0.5 * params.rho * pi * B .* (RO.^4 - params.rs^4); % Inertia grid
    
    % Plot Settings
    figure('Color', 'w');
    
    % Volume Contours
    [C1, h1] = contour(RO*1000, B*1000, Vol_grid, 20, 'LineWidth', 1.2);
    ylabel(colorbar, 'Flywheel Volume (m^3)', 'FontSize', 11, 'FontWeight', 'bold');
    hold on;
    clabel(C1, h1, 'FontSize', 8, 'Color', 'k');
    
    % Plotting Required Inertia Constraint Boundary (J = J_req)
    [C2, h2] = contour(RO*1000, B*1000, J_grid, [params.J_req, params.J_req], ...
                       'm', 'LineWidth', 1.5);
    
    % Plotting Fixed Thickness Constraint Line (b = 1 inch)
    p_bline = yline(b_target*1000, 'm--', 'LineWidth', 1.5);
    
    % Plot optimization path & iteration steps
    p_steps = plot(x_steps(:,1)*1000, x_steps(:,2)*1000, 'ro--', 'MarkerFaceColor', 'r', 'LineWidth', 1.5); % Step points
    p_initial = plot(x0(1)*1000, x0(2)*1000, 'bo', 'MarkerFaceColor', 'b'); %Initial guess
       
    % Mark Optimal Design Point
    p_opt = plot(ro_opt*1000, b_opt*1000, 'o', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
    
    % Formatting
    xlabel('Outer Radius r_o (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Axial Thickness b (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    axis tight;
    
    legend('Volume Contours', 'Boundary (J_{req} = 0.1008 kg\cdotm^2)', ...
            'Boundary (Axial Thickness = 1 inch)', 'Path Iteration Steps', ...
            'Initial Guess', 'Optimal Design Point', ...
           'Location', 'Northeast');
end

%% Volume function
function V = obj_volume(x, params)
% Defining all variables
    r_o = x(1);
    b   = x(2);

    % Calculating volume
    V = pi * (r_o^2 - params.rs^2) * b;
end

%% Nonlinear Constraint Function: Makes sure c <= 0
function [c, ceq] = nonlcon(x, params)
% Defining all variables
    r_o = x(1); % Outer radius [m]
    b   = x(2); % Axial thickness [m]
    
    % Moment of inertia for solid flywheel with center hole
    J = 0.5 * params.rho * pi * b * (r_o^4 - params.rs^4);
    
    % Nonlinear constraints (c <= 0)
    c = params.J_req - J; % J >= J_req
    ceq = [];
end

%% Function to store iteration history
function stop = outfun(x, optimValues, state)
    stop = false; % Keeps solver running
    global x_steps;
    if strcmp(state, 'iter')
        x_steps = [x_steps; x(1:2)]; %Store outer radius and thickness at each step
    end
end