%% Project 1: Air Compressor Flywheel 
% Goal to systematically determine the ideal size of a flywheel when 
% given air compressor design specifications. The goal of the flywheel 
% design process was to provide a coefficient of speed fluctuation no 
% greater than 5%.

% Clear workspace, command window, and any plots in MATLAB
clear;
close all;
clc;

%% Determine the Vertical Position of the Piston-Crank

% Set parameters
AoA = 0.0508;                       % Crank length, [m]
AB = 0.1524;                        % Connecting rod length, [m]

cyl_dia = 0.0635;                   % Cylinder diameter, [m]
cyl_area = (pi/4).*(cyl_dia.^2);    % Cylinder area, [m^2]
h_clear = 0.0222758;                % Head Clearance, [m]
p_atm = 101325;     % Inlet/Atmospheric Pressure [Pa] absolute, 14.7 [PSI]
p_out_g = 344738;                   % Outlet Pressure, [Pa] gage, 50 [PSI]
p_out_a = p_out_g + p_atm;          % Outlet Pressure, [Pa] absolute
comp_speed = 800*pi/30;             % Compressor Speed, [rad/s]

% Set theta as the angle CCW from bottom dead center from 0 to 2*PI
theta = linspace(0,2*pi,1000000);   % Crank Angle, [rad]

% Calculate phi as angle between crank and connecting rod
phi = asin(AoA.*sin(theta)./AB);    % Rod Angle, [rad]

% Calculate vertical displacement by correlating theta and height
y_B = AB*cos(asin(AoA*sin(theta)/AB)) - AoA*cos(theta); % Height [m]

% Parameterize the maximum and minimum heights of the piston
y_max = max(y_B);   % Maximum piston height [m]
y_min = min(y_B);   % Minimum piston height [m]


%% Plot theta vs height of piston
%{
% Plot theta vs piston height
plot(theta, y_B,'k-',linewidth=1)
%title("Displacement vs Crank Angle")
xlabel("Crank Angle, \theta (rad)");
ylabel("Vertical Position (m)");
xlim([0,max(theta)]);
ylim([min(y_B)*0.9,max(y_B)*1.1]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
% Plot centerline (average height)
yline((y_max+y_min)/2,'linestyle','--',Label="Centerline",linewidth=1);
grid on;
hold off;

% Export plot
fig = gcf;
exportgraphics(fig, 'Displacement_vs_Angle.png', 'Resolution', 300);
%}

% Parameterize top dead center and bottom dead center
TDC = y_max;            % Top Dead Center [m]
BDC = y_min;            % Bottom Dead Center [m]

% Print results of max and min piston position
fprintf("Piston height range from [%0.3f m,%0.3f m], or " + ...
    "[%0.1f in, %0.1f in]\n",y_min,y_max, y_min*39.3701, y_max*39.3701)
fprintf("Top Dead Center = %0.2f m\nBottom Dead Center = %0.2f m\n\n" + ...
    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n",TDC,BDC)

%% Thermo analysis, determine pressure as a function of volume
% volume is a function of height which is a function of theta
% Pv^n = const, Polytropic expansion and compression, 1->2 and 3->4
n = 1.3;                                    % Given polytropic constant
vol_height = y_max + h_clear - y_B;         % Height of air in cylinder [m]
vol_theta = cyl_area .* vol_height;         % Volume in cylinder [m^3]
vol_BDC = max(vol_theta);                   % Cylinder volume at BDC [m^3]
vol_TDC = min(vol_theta);                   % Cylinder volume at TDC [m^3]

% Calculate constant compression value of Pv^n
Pvn_constant_compression = p_atm .* (vol_BDC).^n;

% Calculate constant expansion value of Pv^n 
Pvn_constant_expansion = p_out_a .* (vol_TDC).^n;

%% Calculate Pressure and volume as a function of theta for entire process
% 1. Volume at the end of polytropic compression (P reaches P_out_a)
V_2 = vol_BDC * (p_atm / p_out_a)^(1/n);

% Piston height corresponding to V_2
y_2 = y_max + h_clear - (V_2 / cyl_area);

% Solve for theta in compression stroke (between 0 and pi)
% Using 1D numerical interpolation over the first half cycle (0 to pi):
idx_half1 = 1:floor(length(theta)/2);
% Interpolates to find what Y corresponds to X_target interp1(X,Y,X_target)
th_2 = interp1(y_B(idx_half1), theta(idx_half1), y_2);

% 2. Volume at the end of polytropic expansion (P drops back to P_atm)
V_4 = vol_TDC * (p_out_a / p_atm)^(1/n);

% Piston height corresponding to V_4
y_4 = y_max + h_clear - (V_4 / cyl_area);

% Solve for theta in expansion stroke (between pi and 2*pi):
idx_half2 = floor(length(theta)/2)+1 : length(theta);
th_4 = interp1(y_B(idx_half2), theta(idx_half2), y_4);

% Print calculated transition angles
fprintf('Compression ends / Discharge starts at theta = %0.4f rad (%0.2f deg)\n', ...
    th_2, rad2deg(th_2));
fprintf('Expansion ends / Intake starts at theta = %0.4f rad (%0.2f deg)\n\n', ...
    th_4, rad2deg(th_4));

% Pre-allocate pressure matrix to fill with pressure values
P_theta = zeros(size(theta)); 

% for loop
for ii = 1:length(theta)
    % Set current crank angle [rad]
    th = theta(ii); 
    % Calculate current piston height from ground [m]
    y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th);
    vol_height_th = y_max + h_clear - y_B; % Volume height at theta(ii)
    vol_th = cyl_area .* vol_height_th; % Volume at theta(ii)

    % Polytropic compression until P = 50 PSI
    if (th >= 0 && th <= th_2) 
        P_theta(ii) = Pvn_constant_compression ./ (vol_th .^ n); 
        
    % Constant pressure discharge after P = 50 PSI
    elseif (th > th_2 && th <= pi) 
        P_theta(ii) = p_out_a;
        
    % Polytropic expansion
    elseif (th > pi && th <= th_4) 
        P_theta(ii) = Pvn_constant_expansion ./ (vol_th .^ n); 
    
    % Constant pressure intanke
    elseif (th > th_4 && th <= 2*pi) 
        P_theta(ii) = p_atm;
    end
end

% P_theta now contains pressure(theta) across the whole process

%% Plot of pressure as a function of theta
%{
figure();
plot(theta, P_theta,'k-','linewidth',1)
%title({'Pressure in Cylinder vs Angle of Crank';''})
xlabel("Crank Angle, \theta (rad)");
ylabel("Pressure (Pa)");
xlim([0,max(theta)]);
ylim([50000,500000]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
grid on;

% Export plot
fig = gcf;
exportgraphics(fig, 'Pressure_vs_Angle.png', 'Resolution', 300);

%% Plot Pressure vs Volume
figure()
plot(vol_theta, P_theta,'k-','linewidth',1)
%title({'Pressure in Cylinder vs Volume';''})
xlabel("Volume (m^3)");
ylabel("Pressure (Pa)");
xlim([0.00005 0.0004]);
ylim([50000,500000]);
grid on;

% Export plot
fig = gcf;
exportgraphics(fig, 'Pressure_vs_Volume.png', 'Resolution', 300);
%}
%% Calculate loading torque
P_theta_gage=P_theta-p_atm; % Adjust P_theta to be gage pressure
% Calculate Force on piston face at each angle theta [N]
F_theta = P_theta_gage.*cyl_area; 
% Calculate Loading torque [Nm] as a function of theta
T_load_theta = F_theta.*AoA.*(sin(theta)-(tan(phi).*cos(theta))); 

% Index the value corresponding to theta = 2pi
idx_2pi = find(theta <= 2*pi, 1, 'last');

% Calculate average torque over one full crank rotation (0 to 2pi)
T_avg = (1 / (2*pi)) * trapz(theta(1:idx_2pi), T_load_theta(1:idx_2pi));
% Print average torque value
fprintf('Average Torque is %0.3f Nm\n', T_avg);

% Plot torque vs theta
%{
figure()
plot(theta,T_load_theta,'k-','linewidth',1)
%title("Torque vs Crank Angle")
xlabel("Crank Angle, \theta (rad)");
ylabel('Torque (Nm)')
% Plot average torque
yline(T_avg,'label','Average Torque','linestyle','--',linewidth=1)
xlim([0,max(theta)])
grid on;
hold on;
%}
%% Calculate where T_load_theta intersects with T_avg
% Difference vector
T_diff = T_load_theta - T_avg;

% Find where the sign changes (zeros)
cross_idx = find(T_diff(1:end-1) .* T_diff(2:end) <= 0);

% Use index to find what theta it intersects at
theta_intersects = zeros(size(cross_idx));
for k = 1:length(cross_idx)
    theta_intersects(k) = theta(cross_idx(k));
end

% Display intersection points
fprintf('Angle Theta at Intersection:\n%0.4f rad & %0.4f rad\n\n', ...
    theta_intersects);

% Find array indices corresponding to two bounding intersection angles
Theta_1 = cross_idx(1);
Theta_2 = cross_idx(2);

% Plot Theta_1 and Theta_2 on Torque vs Theta graph as bounds
%{
plot(theta(Theta_1), T_avg,'o','markerfacecolor','k','markeredgecolor','k')
plot(theta(Theta_2),T_avg,'o','markerfacecolor','k','markeredgecolor','k')

% Export plot
fig = gcf;
exportgraphics(fig, 'Torque_vs_Angle.png', 'Resolution', 300);
%}
%% Find Kinetic Energy
% Integrate delta torque to find peak kinetic energy fluctuation [J]
dKE = trapz(theta(Theta_1:Theta_2), T_load_theta(Theta_1:Theta_2) - T_avg);

% Display change in KE
fprintf('Maximum Change in Kinetic Energy = %0.3f J\n', dKE);

% Choose cf = 0.05 (maximum allowable) and find related rotational Inertia
cf = 0.05;
J_cf = dKE./(cf.*comp_speed.^2);

%% Flywheel Values
% Parameterize knowns
rho = 7200;             % material density, [kg/m^3]
b = 0.0254;             % axial thickness, [m]
t_t = 0.0381;           % Radial Thickness, [m]

% Calculate Unknowns by varying outer diameter
D_o = linspace(0.1,1,100000);               % Outer Diameter, [m]
D_i = D_o - 2*t_t;                          % Inner Diameter, [m]
r_o = D_o./2;                               % Outer Radius, [m]
r_i = D_i./2;                               % Inner Radius, [m]
mass = rho.*pi.*b.*(r_o.^2-r_i.^2);         % mass, [kg]
J_D = 1/2.*mass.*(r_o.^2+r_i.^2);           % Rotational Inertia, [kg/m^2]


%% Find diameters, rotational inertia, and mass
% Find outer radius when J is J_cf 
% Difference vector
J_diff = J_D - J_cf;

% Find where the sign changes (zeros)
J_cross_idx = find(J_diff(1:end-1) .* J_diff(2:end) <= 0);

% Use index to find what outer radius it intersects at
% Pre-allocate intersect matrix
r_o_intersects = zeros(size(J_cross_idx));
for k = 1:length(J_cross_idx)
    r_o_intersects(k) = r_o(J_cross_idx(k));
end

d_outer = 2*r_o(J_cross_idx(1));    % Outer diameter, [m]
d_inner = 2*r_i(J_cross_idx(1));    % Inner diameter, [m]
mass_final = mass(J_cross_idx(1));  % Mass from rotational inertia, [kg]

% Display intersection points
fprintf('Outer diameter: %0.4f m\n',d_outer);
fprintf('Inner diameter: %0.4f m\n',d_inner);
fprintf('Mass: %0.4f kg\n',mass_final);
fprintf('Rotational Inertia: %0.4f kg*m^2\n\n',J_cf);

%% Find both maximum and minimum angular velocities
w_avg = comp_speed;              % Nominal average speed [rad/s]
w_max = w_avg * (1 + cf / 2);    % Maximum angular velocity [rad/s]
w_min = w_avg * (1 - cf / 2);    % Minimum angular velocity [rad/s]

% Print max and min angular velocities
fprintf('Max angular velocity: %0.2f rad/s\n', w_max);
fprintf('Min angular velocity: %0.2f rad/s\n', w_min);

%% Calculate and Plot Flywheel Speed vs Crank Angle
% Calculate net torque vector
T_net = T_load_theta - T_avg;

% Define ODE handle using interp1
dwdtheta = @(th,w) interp1(theta,T_net,th,'linear','extrap')/(J_cf*w);

w0 = 82.3;          % Initial guess, found so that w_avg = average w

% Use ode45 to find angular velocity at each theta
[theta_ode, w_ode] = ode45(dwdtheta, [0,2*pi], w0);

% Plot angular velocity vs crank position
%{
figure()
plot(theta_ode, w_ode,'k-', 'LineWidth', 1);
% Plot initial angular velocity
yline(w_avg, 'k--', 'Average Speed (\omega_0)', linewidth=1, ...
    LabelHorizontalAlignment='center');
%title('Flywheel Speed vs Crank Angle');
xlabel("Crank Angle, \theta (rad)");
ylabel('Rotational Velocity, \omega (rad/s)');
xlim([0, 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi","4\pi"]);
grid on;

% Export plot
fig = gcf;
exportgraphics(fig, 'Speed_vs_Angle.png', 'Resolution', 300);
%}
%% Find Required Motor Power
% Print power for motor with and without flywheel for comparison
fprintf('\nPower needed for motor without the flywheel = %0.2f W\n', ...
    max(T_load_theta)*w_avg)
fprintf('Power needed for motor with the flywheel = %0.2f W\n', ...
    mean(T_load_theta)*w_avg)

%% Find Work and Efficiency
% Set bounds to constant pressure discharge where theta >th_2 and <=pi
Bounds = (theta > th_2) & (theta <= pi);
P_bounded = P_theta_gage(Bounds);
vol_bounded = vol_theta(Bounds);

% Calculate work out at constant pressure discharge
W_out = trapz(vol_bounded,P_bounded);       % Work out [J]

% Calculate work in for system
Power_in = T_avg*w_avg;     % Power in [W]
time = 2*pi/comp_speed;     % Time for one cycle [s]
W_in = Power_in*time;       % Work in [J]

% Calculate efficiency W_out/W_in
Efficiency = abs(W_out/W_in);

% Print efficiency value as a percentage
fprintf('Efficiency of air compressor = %0.2f percent\n', Efficiency*100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Homework 2
%% Evaluate W
% Function outputting work saved and total compression work for the 2 stage
% compressor
function [Work_saved, Work_total] = Work_calc(P_int,efficiency,n,vol_BDC, ...
    p_out_a)
    % Finds work saved and total compression work of a 2 stage piston
    Temp_1 = 295.15;            % Ambient air at 20 [C] = 295.15 [K]
    P_1 = 101325;               % Atmospheric Pressure [Pa]
    P_2 = p_out_a;              % Maximum pressure at state 2 [Pa]

    Temp_int = Temp_1.*(P_int/P_1).^((n-1)/n);    
    Temp_3 = Temp_int - efficiency.*(Temp_int-Temp_1);
    
    % Work from state 1 at p_atm to P_int:
    W1 = (n/(n-1)).*P_1.*vol_BDC.*(((P_int./P_1).^((n-1)./n))-1);
    
    % Volume at intermediate stage
    V_int = vol_BDC.*(P_1./P_int).^(1/n).*(Temp_3./Temp_int);

    % Work from P_int to state 2 at p_out_a
    W2 = (n/(n-1)).*P_int.*V_int.*(((P_2./P_int).^((n-1)./n))-1);

    % Calculate Work from state 1 to state 2 uninterrupted
    W12 = (n/(n-1)).*P_1.*vol_BDC.*(((P_2./P_1).^((n-1)./n))-1);

    % Calculate total work and work saved to report as outputs
    Work_total = W1+W2;
    Work_saved = W12-Work_total;
end

P_int = linspace(p_atm, p_out_a,100000);

% Efficiency of 1.0:
[Work_saved_1,Work_total_1] = Work_calc(P_int, 1,n,vol_BDC,p_out_a);
idx_W1 = find(Work_total_1 == min(Work_total_1));
Min_W1 = Work_total_1(idx_W1);
Min_P_int1 = P_int(idx_W1);

% Efficiency of 0.6:
[Work_saved_2,Work_total_2] = Work_calc(P_int, 0.6, n,vol_BDC,p_out_a);
idx_W2 = find(Work_total_2 == min(Work_total_2));
Min_W2 = Work_total_2(idx_W2);
Min_P_int2 = P_int(idx_W2);


figure()
plot(P_int,Work_total_1,'linewidth',1,'Color','k')
hold on;
plot(Min_P_int1,Min_W1,'ko','markerfacecolor','k')
hold on;
plot(P_int,Work_total_2,'k--','linewidth',1)
hold on;
plot(Min_P_int2,Min_W2,'ks','markerfacecolor','k')
xlabel('Intermediate Pressure, P_{int} [Pa, abs]')
ylabel('Total Compression Work, W_{comp,t} [J]')
legend('Efficiency \eta = 1.0','Minimum Work for \eta = 1.0','Efficiency \eta = 0.6', ...
    'Minimum Work for \eta = 0.6')
grid on;

% Print results
fprintf(['\nIntermediate Pressure for minimum calculated work of ' ...
    '%0.2f J at 100%% effectiveness is %0.f Pa\n'],Min_W1,Min_P_int1)
fprintf(['Intermediate Pressure for minimum calculated work of ' ...
    '%0.2f J at 60%% effectiveness is %0.f Pa\n'],Min_W2,Min_P_int2)

fig = gcf;
exportgraphics(fig, 'Work_vs_IntermediatePressure.png', 'Resolution', 300);

%% Optimization

% Define function to call only work saved value
function workVal = compWork(P_int,efficiency,n,vol_BDC, p_out_a)
    [~,workVal] = Work_calc(P_int,efficiency,n,vol_BDC, p_out_a);
end



% begin function to minimize work for 2-stage compressor
function min2stagepcompwork(P_int,Work_total_1,Work_total_2)
    
    p_atm = 101325;     % Recall required input/output values
    p_out_g = 344738;                   
    p_out_a = p_out_g + p_atm;
    n = 1.3;
    vol_BDC = 3.923049646188528 .* 10.^(-4);

% a.) Provide the objective function 
    objFun1  = @(P_int) compWork(P_int,1,n,vol_BDC, p_out_a);
    objFun06 = @(P_int) compWork(P_int,0.6,n,vol_BDC, p_out_a);

% b.) List any and all constraints 
    lb = p_atm;       % lower bound for P_int is intake pressure
    ub = p_out_a;     % upper bound for P_int is outlet pressure
    p_ig = 0.25 .* (p_out_a - p_atm) + p_atm; % Pressure initial guess [Pa]
    % W_ig = 60; % Compression work initial guess [J]
    % initGuesses = [p_ig; W_ig];

% c.) Plot optimization process based on desired initial guess

    % Setup History Tracking and Default Options
    history1.x = p_ig; % Array to record solver trajectory
    effVal = 1;

    options = optimoptions('fmincon', ...
        'Display', 'iter', ...      % Print iterations to command window
        'OutputFcn', @saveHistory, ... % Callback to record x at each step
        'StepTolerance', 1e-12, ...
        'FiniteDifferenceStepSize', 1e-3, ...
        'OptimalityTolerance', 1e-10); 

    [p_int_opt1, minWork1] = fmincon(objFun1, p_ig, [], [], [], [], ...
        lb, ub, [], options);


    % Display results in command window
    fprintf('\n--- 100 Percent Efficiency Optimization Results ---\n');
    fprintf('Optimal Intermediate Pressure: %.4f\n', p_int_opt1);
    fprintf('Minimum Work: %.4f\n', minWork1);
 
    % Plot results
    % P_int = linspace(p_atm, p_out_a, 100000);
    
    % X1 = meshgrid(10^5:10^4:4.5*10^5);
    % Z1  = compWork(P_int,  1,n,vol_BDC, p_out_a); % Evaluate objective over grid
    % Z06 = compWork(P_int,0.6,n,vol_BDC, p_out_a); % Evaluate objective over grid
    Work_history1 = compWork(history1.x, 1,n,vol_BDC, p_out_a);


    % 100% Efficiency Optimization Plot:
    
    figure()
    plot(P_int,Work_total_1,'k-','linewidth',1,'displayname', ...
        'Work for \eta = 1.0')
    hold on;
    % Draw optimization path
    plot(history1.x(1,:), Work_history1, 'ro-', 'LineWidth', 2, ...
        'MarkerFaceColor', 'y', 'MarkerSize', 6, 'DisplayName', ...
        'Solver Trajectory');
    
    hold on;

    % Highlight start and end points
    plot(p_ig, 64.4242, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b', ...
        'DisplayName', 'Start (P_{int,0})');
    plot(p_int_opt1, minWork1, 'go', 'MarkerSize', 6, 'MarkerFaceColor',...
        'g', 'DisplayName', 'Optimum (P_{int}^*)');

    % Axis formatting
    xlabel('Intermediate Pressure, P_{int} [Pa, abs]', 'FontSize', 14); 
    ylabel('Total Compression Work, W_{comp,t} [J]', 'FontSize', 14);
    legend('Location', 'northeast');
    axis([1.85E5 2.15E5 64.24 64.44]);

    fig = gcf;
    exportgraphics(fig, 'Optimized_1.png', 'Resolution', 300);



    % 60 Efficiency
     % Setup History Tracking and Default Options
    history06.x = p_ig; % Array to record solver trajectory
    effVal = 0.6;

    % options = optimoptions('fmincon', ...
    %     'Display', 'iter', ...      % Print iterations to command window
    %     'OutputFcn', @saveHistory, ... % Callback to record x at each step
    %     'StepTolerance', 1e-12, ...
    %     'FiniteDifferenceStepSize', 1e-3, ...
    %     'OptimalityTolerance', 1e-10); 

    [p_int_opt06, minWork06] = fmincon(objFun06, p_ig, [], [], [], [], ...
            lb, ub, [], options);

    % Display results in command window
    fprintf('\n--- 60 Percent Efficiency Optimization Results ---\n');
    fprintf('Optimal Intermediate Pressure: %.4f\n', p_int_opt06);
    fprintf('Minimum Work: %.4f\n', minWork06);

    % Plot results
    % P_int = linspace(p_atm, p_out_a, 100000);
    
    % X1 = meshgrid(10^5:10^4:4.5*10^5);
    % Z1  = compWork(P_int,  1,n,vol_BDC, p_out_a); % Evaluate objective over grid
    % Z06 = compWork(P_int,0.6,n,vol_BDC, p_out_a); % Evaluate objective over grid
    Work_history06 = compWork(history06.x, 0.6,n,vol_BDC, p_out_a);


    % 60% Efficiency Optimization Plot:
    % figure('Name', '60% Efficiency', 'NumberTitle', 'off');

    % Draw optimization path
    figure()
    plot(P_int,Work_total_2,'k-','linewidth',1,'displayname', ...
        'Work for \eta = 0.6')
    hold on;
    plot(history06.x(1,:), Work_history06, 'ro-', 'LineWidth', 2, ...
        'MarkerFaceColor', 'y', 'MarkerSize', 6, 'DisplayName', ...
        'Solver Trajectory');
    
    hold on;

    % Highlight start and end points
    plot(p_ig, 66.7522, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b', ...
        'DisplayName', 'Start (P_{int,0})');
    plot(p_int_opt06, minWork06, 'go', 'MarkerSize', 6, 'MarkerFaceColor',...
        'g', 'DisplayName', 'Optimum (P_{int}^*)');

    % Axis formatting
    xlabel('Intermediate Pressure, P_{int} [Pa, abs]', 'FontSize', 14); 
    ylabel('Total Compression Work, W_{comp,t} [J]', 'FontSize', 14);
    legend('Location', 'northeast');
    axis([1.85E5,2.15E5, 66.64,66.76]);

    fig = gcf;
    exportgraphics(fig, 'Optimized_2.png', 'Resolution', 300);




    % Nested Output Function to store steps
    function stop = saveHistory(x, optimValues, state)
        stop = false;
        if strcmp(state, 'iter') || strcmp(state, 'init')
            if effVal == 1
                history1.x = [history1.x, x];
            else
                history06.x = [history06.x, x];
            end
        end
    end
end


min2stagepcompwork(P_int,Work_total_1,Work_total_2)