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
comp_speed = 84;               % Compressor Speed, [rad/s]

theta=linspace(0,2*pi,10000);       % Crank Angle, [rad]
phi = asin(AoA.*sin(theta)./AB);   % Rod Angle, [rad]

% Calculate theta to height correlation
y_B = AB*cos(asin(AoA*sin(theta)/AB)) - AoA*cos(theta);

y_max = max(y_B);
y_min = min(y_B);


%% Plot theta vs height of piston
subplot(1,2,1)
plot(theta, y_B,linewidth=1)
title("Height of Piston vs Angle of Crank")
xlabel("\theta (rad)");
ylabel("y position (m)");
xlim([0,max(theta)]);
ylim([min(y_B)*0.9,max(y_B)*1.1]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
yline((y_max+y_min)/2,'linestyle','--',Label="Centerline");
grid on;
hold off;

%% Plot phi vs height of piston
subplot(1,2,2)
plot(phi, y_B,linewidth=1)
title("Height of Piston vs Angle of Rod")
xlabel("\phi (rad)");
ylabel("y position (m)");
xlim([min(phi)*1.1,max(phi)*1.1]);
ylim([min(y_B)*0.9,max(y_B)*1.1]);
xticks(0:pi/2:max(phi))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
grid on;


% Calculate top dead center and bottom dead center
TDC = y_max;          % Top Dead Center
BDC = y_min;  % Bottom Dead Center

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
    th = theta(ii); %current crank angle [rad]
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

   elseif (th >= 2*pi && th <= 2.39 + 2*pi) %polytropic compression
      P_theta(ii) = Pvn_constant_compression ./ (vol_th .^ n);

   elseif (th > 2.39 + 2*pi && th <= 3*pi) %constant pressure discharge
      P_theta(ii) = p_out_a;

   elseif (th > 3*pi && th <= 4.48 + 2*pi) %polytropic expansion
      P_theta(ii) = Pvn_constant_expansion ./ (vol_th .^ n); 

   else %constant pressure intake
      P_theta(ii) = p_atm;

   end
end

%% Plot of Pressure as a function of Theta
figure;
subplot(1,2,1)
plot(theta, P_theta,'linewidth',1)
title({'Pressure in Cylinder vs Angle of Crank';''})
xlabel("\theta (rad)");
ylabel("Pressure (Pa)");
xlim([0,max(theta)]);
ylim([50000,500000]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
grid on;

%% Plot Pressure vs Volume
subplot(1,2,2)
plot(vol_theta, P_theta,'linewidth',1)
title({'Pressure in Cylinder vs Volume';''})
xlabel("Volume (m^3)");
ylabel("Pressure (Pa)");
xlim([0.00005 0.0004]);
ylim([50000,500000]);
grid on;


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
plot(theta,T_load_theta,'linewidth',1)
title("Torque vs. Theta")
xlabel('\theta (rad)')
ylabel('Torque (Nm)')
yline(T_avg,'label','Average Torque','linestyle','--')
xlim([0,max(theta)])
grid on;


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
w_avg = 800;
J_guess = dKE./(cf_guess.*w_avg.^2);

% Range J around estimate and find resulting ranged cf
J_range = linspace(0,J_guess*2,10000);
cf = dKE./(w_avg.^2.*J_range);

% Find where cf is minimized and print minimum cf
[min_cf, idx] = min(abs(cf));
fprintf('The minimum cf is %0.4f\n',min_cf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Flywheel
%{
rho = 7200;             % material density, [kg/m^3]
b = 0.0254;             % axial thickness, [m]
t_t = 0.0381;           % Radial Thickness, [m]

% Calculate Unknowns
D_o = linspace(0.1,0.25);                   % Outer Diameter, [m]
D_i = D_o - 2*t_t;                          % Inner Diameter, [m]
r_o = D_o./2;                               % Outer Radius, [m]
r_i = D_i./2;                               % Inner Radius, [m]
mass = rho.*pi.*t_t.*(r_o.^2-r_i.^2);       % mass, [kg]
J = 1/2.*mass.*(r_o.^2-r_i.^2);             % Rotational Inertia, [kg/m^2]
%}

%% Final Project Goal, 'Coefficient Fluctuation in Speed' C_f
% we should already have values of omega at this point

% w_max = max(w);
% w_min = min(w); 
% w_avg = mean(w);

%C_f = (w_max - w_min) ./ w_avg;

%if C_f <= 0.05
   % disp('C_f value is sufficient')

%else
  %  disp('Warning: C_f value is too large')

%end

%% Calculate flywheel inertia
%Use both inertia equations to relate to flywheel specifications
%J = (dKE / )

