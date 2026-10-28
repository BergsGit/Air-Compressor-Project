%% Homework 2, Problem 2, Part 2
% Goal: to optimize the design of a spoked flywheel to minimize mass
% while also meeting design requirements and parameters.

% Clear workspace, command window, and any plots in MATLAB
clc
clear
close all 

%% Setup Parameters
target_J = 0.1008; % Required J [kg*m^2]
b = 0.0254; % Axial thickness [m] (1 inch)
t_t = 0.0381; % Radial thickness [m] (1.5 inches)

% Define constants in a structure
params.rho = 7150;                % Material density [kg/m^3]
params.rs = 0.012;                % Shaft radius [m] (12 mm)
params.rhub = 0.034;              % Hub outer radius [m] (34 mm)
params.J_req = target_J;          % Required moment of inertia [kg*m^2]
params.b = b;                     % Axial thickness [m]
params.t_t = t_t;                 % Radial thickness [m]

%% Solve solid flywheel
% Required outer radius value
ro_solid = (((2 * target_J)/(params.rho * pi * b)) + (params.rs^4))^(1/4); 
% Required mass value
mass_solid = (2 * target_J)/((ro_solid^2) + (params.rs^2));

% Call the optimization function
flywheel_optimization(params);

%% Optimization Function Definition
function flywheel_optimization(params)
    % Unpack structure
    rho     = params.rho;
    rs      = params.rs;
    rhub    = params.rhub;
    target_J = params.J_req;
    t_rad   = params.t_t;
    b_axial = params.b;
    
    % Design Constraints
    max_ro  = 0.2032;   % Max outer radius (8 inches)
    min_Ws  = 0.010;    % Min feature size (10 mm)
    
    % Global array to store iteration history
    global x_steps;
    
    % Initial guess for fmincon
    x0 = [0.1, 0.015, b_axial, t_rad]; % Initial guess [r_o, Ws, b_axial, t_rad]
    x_steps = x0(1:2); % Initialize itertation history with starting point
    
    % Defining upper and lower bounds
    lowerbound = [rhub + t_rad, min_Ws, b_axial, t_rad];
    upperbound = [max_ro, 0.07, b_axial, t_rad];
    
    % Optimization Options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'OutputFcn', @outfun);
    
    % Run fmincon
    [x_opt, V_opt, exitflag, output] = fmincon(@(x) obj_volume(x, params), x0, [], [], [], [], ...
                             lowerbound, upperbound, @(x) nonlcon(x, params), options);
    
    % Store optimized results
    ro_opt    = x_opt(1); % Optimized outer radius [m]
    Ws_opt    = x_opt(2); % Optimized spoke width [m]
    b_opt     = x_opt(3); % Optimized axial thickness [m]
    t_rad_opt = x_opt(4); % Optimized radial thickness [m]
    ri_opt    = ro_opt - t_rad_opt; % Optimized inner radius [m]
    mass_opt  = V_opt * rho; % Optimized mass [kg]
    
    % Print Resultss
    fprintf('Optimal Spoked Flywheel Design Results\n');
    fprintf('Iterations: %d\n', output.iterations);
    fprintf('Function Evaluations: %d\n', output.funcCount);
    fprintf('Total Mass: %.2f kg\n', mass_opt);
    
    %% Generate Contour Plot over (r_o, Ws) Space
    ro_grid = linspace(rhub + t_rad, max_ro, 200); % Outer Radius values
    Ws_grid = linspace(0, 0.07, 200);  % Spoke width values
    [RO, WS] = meshgrid(ro_grid, Ws_grid);
   
    % Evaluate Volume and Inertia across the domain grid
    V_grid = zeros(size(RO)); % Empty volume grid
    t_req_grid = zeros(size(RO)); % Empty thickness grid
    
    % Loop to calculate properties for entire grid
    for i = 1:size(RO, 1) 
        for j = 1:size(RO, 2)
            % Variables for loop run (i)
            ro_value = RO(i,j); % Outer radius value (i)
            Ws_value = WS(i,j); % Spoke width value (i)
            ri_value = ro_value - t_rad; % Inner radius value (i)
            
            % Calculate properties per unit thickness using properties function
            [V_prime, J_prime] = properties(ro_value, ri_value, Ws_value, params);
            
            % Required thickness to achieve required J
            t_req = params.J_req / J_prime;
            t_req_grid(i,j) = t_req;
            V_grid(i,j) = V_prime * t_req;
        end
    end
    
    % Plot Settings
    figure('Color', 'w');
    
    % Volume Contours
    [C1, h1] = contour(RO*1000, WS*1000, V_grid*1e6, 15, 'LineWidth', 1.2);
    ylabel(colorbar, 'Total Volume (cm^3)', 'FontSize', 11, 'FontWeight', 'bold');
    hold on;
    
    % Plot optimization path & iteration steps
    p_steps = plot(x_steps(:,1)*1000, x_steps(:,2)*1000, 'ro--', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'LineWidth', 2); % Step points
    p_initial = plot(x0(1)*1000, x0(2)*1000, 'bo', 'MarkerFaceColor', 'b'); % Initial guess
  
    % Label values on contour lines
    clabel(C1, h1, 'FontSize', 8, 'Color', 'k');
    
    % Constraint Boundaries
    p_bnd1 = contour(RO*1000, WS*1000, t_req_grid*1000, [25.4 25.4], 'm--', 'LineWidth', 1.5);
    p_bnd2 = line([min(ro_grid)*1000, max(ro_grid)*1000], [10 10], 'Color', 'm', 'LineWidth', 1.5);
    
    % Mark Optimal Design Point
    p_opt = plot(ro_opt*1000, Ws_opt*1000, 'o', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
    
    % Formatting
    xlabel('Outer Radius r_o (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Spoke Width W_s (mm)', 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    axis tight;

    legend('Volume Contours', 'Path Iteration Steps', ...
             'Initial Guess', 'Boundary (t = 1 in)', ...
             'Boundary( W_s \geq 10 mm)', 'Optimal Design Point', ...
            'Location', 'Northeast');
end

%% Volume function
function V = obj_volume(x, params)
% Defining all variables
    r_o = x(1); % Outer radius [m]
    Ws = x(2); % Spoke width [m]
    b = x(3); % Axial thickness [m]
    t_rad = x(4); % Radial thickness [m]
    r_i = r_o - t_rad; %Inner radius [m]
   
    % Calculating component volumes 
    V_hub = pi * (params.rhub^2 - params.rs^2) * b; % Volume of hub
    V_rim = pi * (r_o^2 - r_i^2) * b; % Volume of rim
    V_spokes = 4 * Ws * (r_i - params.rhub) * b; % Volume of spokes
    
    V = V_hub + V_rim + V_spokes; % Total volume
end

%% Nonlinear Constraint Function: Makes sure c <= 0
function [c, ceq] = nonlcon(x, params)
% Defining all variables
    r_o = x(1); % Outer radius [m]
    Ws = x(2); % Spoke width [m]
    b = x(3); % Axial thickness [m]
    t_rad = x(4); % Outer radius [m]
    ri = r_o - t_rad; %Inner radius [m]
    
    % Unpack parameters
    rho   = params.rho;
    rs    = params.rs;
    rhub  = params.rhub;
    J_req = params.J_req;
    
    % Inertia calculations for each component of flywheel
    J_hub = 0.5 * rho * pi * b * (rhub^4 - rs^4); % Hub inertia
    J_rim = 0.5 * rho * pi * b * (r_o^4 - ri^4); % Rim inertia
    m_spoke = rho * Ws * b * (ri - rhub); % Mass of one spoke
    J_spoke = m_spoke * ((1/3)*(ri^2 + rhub*ri + rhub^2) + (1/12)*Ws^2); % Single spoke inertia
    
    J_total = J_hub + J_rim + 4 * J_spoke; % Overall inertia
    
    % Nonlinear constraints (c <= 0)
    c(1) = J_req - J_total; % J >= J_req
    c(2) = rhub - ri; % r_i >= r_h
    ceq = [];
end

%% Properties per Thickness Function
% Computes volume and moment of inertia per unit thickness (1 meter)
function [V_prime, J_prime] = properties(ro, ri, Ws, params)
% Unpack parameters
    rho  = params.rho;
    rs   = params.rs;
    rhub = params.rhub;
    
    % X-section volume per unit thickness calculation
    V_prime = pi * (rhub^2 - rs^2) + pi * (ro^2 - ri^2) + 4 * Ws * (ri - rhub);
    
    % Inertia per unit thickness calculation
    J_hub_prime = 0.5 * rho * pi * (rhub^4 - rs^4); % Hub inertia/m
    J_rim_prime = 0.5 * rho * pi * (ro^4 - ri^4); % Rim inertia/m
    m_spoke_prime = rho * Ws * (ri - rhub); % Mass of one spoke/m
    J_spoke_prime = m_spoke_prime * ((1/3)*(ri^2 + rhub*ri + rhub^2) + (1/12)*Ws^2); %Inertia of one spoke/m
    
    J_prime = J_hub_prime + J_rim_prime + 4 * J_spoke_prime; % Overall inertia/m
end

%% Function to store iteration history
function stop = outfun(x, optimValues, state)
    stop = false; % Keeps solver running
    global x_steps;
    if strcmp(state, 'iter')
        x_steps = [x_steps; x(1:2)]; % Store outer radius and spoke width at each step
    end
end