clear;
close all;
clc;

%% Position of a Slider-Crank

AoA = 0.0508;                  % Crank length, [m]
AB = 0.1524;                   % Connecting rod length, [m]

cyl_dia = 0.0635;              % Cylinder diameter, [m]
cyl_area = (pi/4).*(cyl_dia.^2);    % Cylinder area, [m^2]
h_clear = 0.0222758;           % Head Clearance, [m]
p_atm = 101325;                % Inlet/Atmospheric Pressure [Pa] absolute, 14.7 [PSI]
p_out_g = 344738;              % Outlet Pressure, [Pa] gage, 50 [PSI]
p_out_a = p_out_g + p_atm;     % Outlet Pressure, [Pa] absolute
comp_speed = 800*pi/30;               % Compressor Speed, [rad/s]

theta=linspace(0,2*pi,1000000);       % Crank Angle, [rad]
phi = asin(AoA.*sin(theta)./AB);   % Rod Angle, [rad]

% Calculate theta to height correlation
y_B = AB*cos(asin(AoA*sin(theta)/AB)) - AoA*cos(theta);

y_max = max(y_B);
y_min = min(y_B);


%% Plot theta vs height of piston
plot(theta, y_B,'k-',linewidth=1)
title("Displacement vs Crank Angle")
xlabel("\theta (rad)");
ylabel("y position (m)");
xlim([0,max(theta)]);
ylim([min(y_B)*0.9,max(y_B)*1.1]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
yline((y_max+y_min)/2,'linestyle','--',Label="Centerline",linewidth=1);
grid on;
hold off;

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


% Calculate top dead center and bottom dead center
TDC = y_max;          % Top Dead Center
BDC = y_min;            % Bottom Dead Center

%% Print results of max and min piston position
fprintf("Piston height range from [%0.3f m,%0.3f m], or " + ...
    "[%0.1f in, %0.1f in]\n",y_min,y_max, y_min*39.3701, y_max*39.3701)
fprintf("Top Dead Center = %0.2f m\nBottom Dead Center = %0.2f m\n\n" + ...
    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n",TDC,BDC)

%% Thermo analysis, determine pressure as a function of volume
% volume is a function of height which is a function of theta
% Pv^n = const, Polytropic expansion and compression, 1->2 and 3->4
n = 1.3;                                    % Given polytropic constant
vol_height = y_max + h_clear - y_B;         % Height of air's volume in the cylinder [m]
vol_theta = cyl_area .* vol_height;         % Volume in cylinder as a function of theta [m^3]
vol_BDC = max(vol_theta);                   % Cylinder volume at BDC [m^3]
vol_TDC = min(vol_theta);                   % Cylinder volume at TDC [m^3]
Pvn_constant_compression = p_atm .* (vol_BDC).^n;% constant compression value of Pv^n
Pvn_constant_expansion = p_out_a .* (vol_TDC).^n;% constant expansion value of Pv^n 

P_theta = zeros(size(theta)); %empty matrix to fill with pressure values

%% loop calculating Pressure and volume as a function of theta for entire thermo process
for ii = 1:length(theta)
    th = theta(ii); % current crank angle [rad]
    y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th);   % Current piston height from gnd [m]
    vol_height_th = y_max + h_clear - y_B; 
    vol_th = cyl_area .* vol_height_th;  

   if (th >= 0 && th <= 2.39) %polytropic compression until P = 50PSI
        P_theta(ii) = Pvn_constant_compression ./ (vol_th .^ n); 
        
   elseif (th > 2.39 && th <= pi) %constant pressure discharge after P = 50 PSI
        P_theta(ii) = p_out_a;
        
   elseif (th > pi && th <= 4.48) %polytropic expansion
        P_theta(ii) = Pvn_constant_expansion ./ (vol_th .^ n); 
    
   elseif (th > 4.48 && th <= 2*pi) %constant pressure intanke
        P_theta(ii) = p_atm;
   end
end


%% Plot of Pressure as a function of Theta
figure();
plot(theta, P_theta,'k-','linewidth',1)
title({'Pressure in Cylinder vs Angle of Crank';''})
xlabel("\theta (rad)");
ylabel("Pressure (Pa)");
xlim([0,max(theta)]);
ylim([50000,500000]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
grid on;

fig = gcf;
exportgraphics(fig, 'Pressure_vs_Angle.png', 'Resolution', 300);

%% Plot Pressure vs Volume
figure()
plot(vol_theta, P_theta,'k-','linewidth',1)
title({'Pressure in Cylinder vs Volume';''})
xlabel("Volume (m^3)");
ylabel("Pressure (Pa)");
xlim([0.00005 0.0004]);
ylim([50000,500000]);
grid on;

fig = gcf;
exportgraphics(fig, 'Pressure_vs_Volume.png', 'Resolution', 300);

%% Calculate loading torque
% Calculate F from thermo equations
F_theta = P_theta.*cyl_area; % Force on piston face at each angle theta [N]
T_load_theta = F_theta.*AoA.*(sin(theta)-(tan(phi).*cos(theta))); % Loading torque, [Nm], f(theta)

% Find the values corresponding to theta = 0 and theta = 2pi
idx_2pi = find(theta <= 2*pi, 1, 'last');

% Calculate average torque over one full crank rotation (0 to 2pi)
T_avg = (1 / (2*pi)) * trapz(theta(1:idx_2pi), T_load_theta(1:idx_2pi));
fprintf('Average Torque is %0.3f Nm\n', T_avg);

figure()
plot(theta,T_load_theta,'k-','linewidth',1)
title("Torque vs Crank Angle")
xlabel('\theta (rad)')
ylabel('Torque (Nm)')
yline(T_avg,'label','Average Torque','linestyle','--',linewidth=1)
xlim([0,max(theta)])
grid on;

fig = gcf;
exportgraphics(fig, 'Torque_vs_Angle.png', 'Resolution', 300);


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
fprintf('Angle Theta at Intersection:\n%0.4f rad & %0.4f rad\n\n',theta_intersects);

% Find array indices corresponding to two bounding intersection angles
Theta_1 = cross_idx(1);
Theta_2 = cross_idx(2);

% Integrate delta torque to find peak kinetic energy fluctuation [J]
dKE = trapz(theta(Theta_1:Theta_2), T_load_theta(Theta_1:Theta_2) - T_avg);
fprintf('Maximum Change in Kinetic Energy = %0.3f J\n', dKE);

% Guess cf to find good estimate for J
cf_guess = 0.05;
J_guess = dKE./(cf_guess.*comp_speed.^2);

% Range J around estimate and find resulting ranged cf
J_range = linspace(0,J_guess,1000000);
cf = dKE./(comp_speed.^2.*J_range);

% Find where cf is minimized and print minimum cf anf resulting J
[min_cf, idx] = min(abs(cf));
J_cf_min = J_range(idx);
fprintf('The minimum cf is %0.2f at J = %0.4f kg*m^2\n\n', min_cf, J_cf_min)

% Find resulting radius and diameter from J
%% Flywheel
rho = 7200;             % material density, [kg/m^3]
b = 0.0254;             % axial thickness, [m]
t_t = 0.0381;           % Radial Thickness, [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flywheel
%(b*Do^3)/2+((b^2)*Do^2)/2+b^4=ro4_ri4;
% Calculate Unknowns
D_o = linspace(0.1,1,100000);               % Outer Diameter, [m]
D_i = D_o - 2*t_t;                          % Inner Diameter, [m]
r_o = D_o./2;                               % Outer Radius, [m]
r_i = D_i./2;                               % Inner Radius, [m]
mass = rho.*pi.*b.*(r_o.^2-r_i.^2);         % mass, [kg]
J = 1/2.*mass.*(r_o.^2+r_i.^2);             % Rotational Inertia, [kg/m^2]

% % Plot J vs outer radius
% figure()
% plot(r_o, J,'linewidth',1)
% title('Rotational Inertia vs Outer Radius')
% yline(J_cf_min,'label','J_{cf min}')
% xlabel('r_o [m]')
% ylabel('J [kg/m^2]')
% grid on;

% Find outer radius when J is J_cf_min 
% Difference vector
J_diff = J - J_cf_min;

% Find where the sign changes (zeros)
J_cross_idx = find(J_diff(1:end-1) .* J_diff(2:end) <= 0);

% Use index to find what ro it intersects at
r_o_intersects = zeros(size(J_cross_idx));
for k = 1:length(J_cross_idx)
    r_o_intersects(k) = r_o(J_cross_idx(k));
end

d_outer = 2*r_o(J_cross_idx(1));
d_inner = 2*r_i(J_cross_idx(1));
mass_final = mass(J_cross_idx(1));

% Display intersection points
fprintf('Outer diameter: %0.4f m\n',d_outer);
fprintf('Inner diameter: %0.4f m\n',d_inner);
fprintf('Mass: %0.4f kg\n',mass_final);
fprintf('Rotational Inertia: %0.4f kg*m^2\n\n',J_cf_min);

% Find both maximum and minimum angular velocities
w_avg = comp_speed;                  % Nominal average speed [rad/s]
w_max = w_avg * (1 + min_cf / 2);    % Maximum angular velocity [rad/s]
w_min = w_avg * (1 - min_cf / 2);    % Minimum angular velocity [rad/s]

% Verify coefficient of fluctuation
C_f_actual = (w_max - w_min) / w_avg;
fprintf('Max angular velocity: %0.2f rad/s\n', w_max);
fprintf('Min angular velocity: %0.2f rad/s\n', w_min);

P_avg = T_avg*w_avg;

%% Calculate and Plot Flywheel Speed vs Crank Angle
% Calculate net torque vector
T_net = T_load_theta - T_avg;

% Define ODE handle using interp1
% 'linear' interpolation with 'extrap' ensures it handles boundary floating-point edge cases
J_flywheel = J_cf_min;
dwdtheta = @(th, w) interp1(theta, T_net, th, 'linear', 'extrap') / (J_flywheel * w);

% Set initial condition and integrate
w0 = 82.3;                    % Initial guess, set to yield w_avg = average w
theta_span = [0, 2*pi];       % Integration span

[theta_ode, w_ode] = ode45(dwdtheta, theta_span, w0);

% Plot w vs theta
figure()
plot(theta_ode, w_ode,'k-', 'LineWidth', 1);
yline(w0, 'k--', 'Average Speed (\omega_0)', linewidth=1, ...
    LabelHorizontalAlignment='center');
title('Flywheel Speed vs Crank Angle');
xlabel('\theta (rad)');
ylabel('\omega (rad/s)');
xlim([0, 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi","7/2\pi","4\pi"]);
grid on;

fig = gcf;
exportgraphics(fig, 'Speed_vs_Angle.png', 'Resolution', 300);

fprintf('\nPower needed for motor without the flywheel = %0.2f W\n', max(T_load_theta)*w_avg)
fprintf('Power needed for motor with the flywheel = %0.2f W\n', mean(T_load_theta)*w_avg)

%% Finding Work and Efficiency

% Const P discharge th > 2.39 && th <= pi
Bounds = (theta > 2.39) & (theta <= pi);
P_bounded = P_theta(Bounds);
vol_bounded = vol_theta(Bounds);
W_out = trapz(vol_bounded,P_bounded);

Power_in = T_avg*w_avg;
time = 2*pi/comp_speed;

W_in = Power_in*time;

Efficiency = abs(W_out/W_in);

fprintf('Efficiency of air compressor = %0.2f percent\n', Efficiency*100)