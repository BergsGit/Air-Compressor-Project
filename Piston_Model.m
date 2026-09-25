clear;
close all;
clc;

%% Position of a Slider-Crank

AoA = 0.0508;           % Crank length, [m]
AB = 0.1524;            % Connecting rod length, [m]

cyl_dia = 0.0635;       % Cylinder diameter, [m]
h_clear = 0.0222758;    % Head Clearance, [m]
p_out = 344738;         % Outlet Pressure, [Pa] gage
comp_speed = 84;        % Compressor Speed, [rad/s]

theta=linspace(0,4*pi,1000);       % Crank Angle, [rad]
phi = asin(AoA.*sin(theta)./AB);   % Rod Angle, [rad]

% Calculate theta to height correlation
y_B = AB*cos(asin(AoA*sin(theta)/AB)) + AoA*cos(theta);

y_max = max(y_B);
y_min = min(y_B);

%% Plot theta vs height of piston
subplot(1,2,1)
plot(theta, y_B,linewidth=1)
title("Height of Piston vs Angle of Crank")
xlabel("\theta (rad)");
ylabel("y position (m)");
xlim([0,max(theta)]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])
yline((y_max+y_min)/2,Label="Centerline");
hold off;

%% Plot phi vs height of piston
subplot(1,2,2)
plot(phi, y_B,linewidth=1)
title("Height of Piston vs Angle of Rod")
xlabel("\phi (rad)");
ylabel("y position (m)");
xlim([min(phi)*1.1,max(phi)*1.1]);
xticks(0:pi/2:max(phi))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])

% Calculate top dead center and bottom dead center
TDC = h_clear;                  % Top Dead Center
BDC = h_clear + (y_max-y_min);  % Bottom Dead Center

%% Print results of max and min piston position
fprintf("Piston height range from [%0.3f m,%0.3f m], or " + ...
    "[%0.1f in, %0.1f in]\n",y_min,y_max, y_min*39.3701, y_max*39.3701)
fprintf("Top Dead Center = %0.2f m\nBottom Dead Center = %0.2f m\n",TDC,BDC)


%% Flywheel
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

%% Calculate loading torque
% Calculate F from thermo equations
% T_load = F*AoA(sin(theta)-tan(phi)cos(theta)); % Loading torque, Nm

