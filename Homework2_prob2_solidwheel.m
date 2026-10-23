clc
clear

%% Homework 2 Solid Design
target_J = 0.1008;
b = 0.0254;
t_t = 0.0381;
% Define constants in a params structure for clean function passing
params.rho = 7150;                % Material density [kg/m^3]
params.rs = 0.012;                % Shaft radius [m] (12 mm)
params.rhub = 0.034;              % Hub outer radius [m] (34 mm)
params.J_req = target_J;          % Required moment of inertia [kg*m^2]
params.b = b;                     % Axial thickness [m]
params.t_t = t_t;                 % Radial thickness [m]
% Solve solid flywheel baseline
ro_solid = (((2 * target_J)/(params.rho * pi * b)) + (params.rs^4))^(1/4); 
mass_solid = (2 * target_J)/((ro_solid^2) + (params.rs^2));
% Display results
fprintf('\nOPTIMAL SOLID FLYWHEEL DESIGN');
fprintf('\nOuter Radius: %.2f mm', ro_solid*1000);
fprintf(' \nTotal Mass: %.2f kg\n', mass_solid);