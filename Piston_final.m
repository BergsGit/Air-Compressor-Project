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

% Plot theta vs piston height
plot(theta, y_B,'k-',linewidth=1)
%title("Displacement vs Crank Angle")
xlabel("\theta (rad)");
ylabel("y position (m)");
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

% %% Plot phi vs height of piston
% subplot(1,2,2)
% plot(phi, y_B,linewidth=1)
% title("Height of Piston vs Angle of Rod")
% xlabel("\phi (rad)");
% ylabel("y position (m)");
% xlim([min(phi)*1.1,max(phi)*1.1]);
% ylim([min(y_B)*0.9,max(y_B)*1.1]);
% xticks(0:pi/2:max(phi))
% xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
%     "7/2\pi", "4\pi"])
% grid on;


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
    if (th >= 0 && th <= 2.39) 
        P_theta(ii) = Pvn_constant_compression ./ (vol_th .^ n); 
        
    % Constant pressure discharge after P = 50 PSI
    elseif (th > 2.39 && th <= pi) 
        P_theta(ii) = p_out_a;
        
    % Polytropic expansion
    elseif (th > pi && th <= 4.48) 
        P_theta(ii) = Pvn_constant_expansion ./ (vol_th .^ n); 
    
    % Constant pressure intanke
    elseif (th > 4.48 && th <= 2*pi) 
        P_theta(ii) = p_atm;
    end
end

% P_theta now contains pressure(theta) across the whole process

%% Plot of pressure as a function of theta
figure();
plot(theta, P_theta,'k-','linewidth',1)
%title({'Pressure in Cylinder vs Angle of Crank';''})
xlabel("\theta (rad)");
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

%% Calculate loading torque
% Calculate Force on piston face at each angle theta [N]
F_theta = P_theta.*cyl_area; 
% Calculate Loading torque [Nm] as a function of theta
T_load_theta = F_theta.*AoA.*(sin(theta)-(tan(phi).*cos(theta))); 

% Index the value corresponding to theta = 2pi
idx_2pi = find(theta <= 2*pi, 1, 'last');

% Calculate average torque over one full crank rotation (0 to 2pi)
T_avg = (1 / (2*pi)) * trapz(theta(1:idx_2pi), T_load_theta(1:idx_2pi));
% Print average torque value
fprintf('Average Torque is %0.3f Nm\n', T_avg);

% Plot torque vs theta
figure()
plot(theta,T_load_theta,'k-','linewidth',1)
%title("Torque vs Crank Angle")
xlabel('\theta (rad)')
ylabel('Torque (Nm)')
% Plot average torque
yline(T_avg,'label','Average Torque','linestyle','--',linewidth=1)
xlim([0,max(theta)])
grid on;
hold on;

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
plot(theta(Theta_1), T_avg,'o','markerfacecolor','k','markeredgecolor','k')
plot(theta(Theta_2),T_avg,'o','markerfacecolor','k','markeredgecolor','k')

% Export plot
fig = gcf;
exportgraphics(fig, 'Torque_vs_Angle.png', 'Resolution', 300);

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

% UNNECESSARY:
% % Plot J vs outer radius
% figure()
% plot(r_o, J,'linewidth',1)
% title('Rotational Inertia vs Outer Radius')
% yline(cf,'label','J_{cf min}')
% xlabel('r_o [m]')
% ylabel('J [kg/m^2]')
% grid on;

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
figure()
plot(theta_ode, w_ode,'k-', 'LineWidth', 1);
% Plot initial angular velocity
yline(w_avg, 'k--', 'Average Speed (\omega_0)', linewidth=1, ...
    LabelHorizontalAlignment='center');
%title('Flywheel Speed vs Crank Angle');
xlabel('\theta (rad)');
ylabel('\omega (rad/s)');
xlim([0, 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi","4\pi"]);
grid on;

% Export plot
fig = gcf;
exportgraphics(fig, 'Speed_vs_Angle.png', 'Resolution', 300);

%% Find Required Motor Power
% Print power for motor with and without flywheel for comparison
fprintf('\nPower needed for motor without the flywheel = %0.2f W\n', ...
    max(T_load_theta)*w_avg)
fprintf('Power needed for motor with the flywheel = %0.2f W\n', ...
    mean(T_load_theta)*w_avg)

%% Find Work and Efficiency
% Set bounds to constant pressure discharge where theta >2.39 and <=pi
Bounds = (theta > 2.39) & (theta <= pi);
P_bounded = P_theta(Bounds);
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