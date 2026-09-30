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

theta=linspace(0,4*pi,1000);       % Crank Angle, [rad]
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

%% Thermo analysis, determine pressure as a function of volume
% volume is a function of height which is a function of theta
% Pv^n = const, Polytropic expansion and compression, 1->2 and 3
n = 1.3;
vol_height = y_max + h_clear - y_B;         % Height of air's volume in the cylinder [m]
vol_total = cyl_area .* vol_height;         % Total volume in the cylinder [m^3]
vol_BDC = max(vol_total);                   % Cylinder volume at BDC [m^3]
vol_TDC = min(vol_total);                   % Cylinder volume at TDC [m^3]
vol_theta = vol_height * cyl_area;          % Volume in cylinder as a fucntion of theta
Pvn_constant_compression = p_atm .* (vol_BDC).^n;% constant compression value of Pv^n
Pvn_constant_expansion = p_out_a .* (vol_TDC).^n;% constant expansion value of Pv^n

P_theta = zeros(size(theta)); %empty matrix to fill with pressure values

%% loop calculating Pressure and volume as a function of theta for entire thermo process
for ii = 1:length(theta)
    th = theta(ii); %current crank angle

   if (th >= 0 && th <= 2.39) %polytropic compression
      y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th); 
      vol_height = y_max + h_clear - y_B; 
      vol_total = cyl_area .* vol_height;  
      P_theta(ii) = Pvn_constant_compression ./ (vol_total .^ n); 
   elseif (th > 2.39 && th <= pi) %constant pressure discharge
        P_theta(ii) = p_out_a;
   elseif (th > pi && th <= 4.48) %polytropic expansion
      y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th);
      vol_height = y_max + h_clear - y_B; 
      vol_total = cyl_area .* vol_height;  
      P_theta(ii) = Pvn_constant_expansion ./ (vol_total .^ n); 
   elseif (th > 4.48 && th <= 2*pi) %constant pressure intanke
       P_theta(ii) = p_atm;
   elseif (th >= 2*pi && th <= 2.39 + 2*pi) %polytropic compression
      y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th);
      vol_height = y_max + h_clear - y_B; 
      vol_total = cyl_area .* vol_height;  
      P_theta(ii) = Pvn_constant_compression ./ (vol_total .^ n);
   elseif (th > 2.39 + 2*pi && th <= 3*pi) %constant pressure discharge
        P_theta(ii) = p_out_a;
   elseif (th > 3*pi && th <= 4.48 + 2*pi) %polytropic expansion
      y_B = AB*cos(asin(AoA*sin(th)/AB)) - AoA*cos(th);
      vol_height = y_max + h_clear - y_B; 
      vol_total = cyl_area .* vol_height;  
      P_theta(ii) = Pvn_constant_expansion ./ (vol_total .^ n); 
   else %constant pressure intake
       P_theta(ii) = p_atm;
   end
end

%% Plot of Pressure as a function of Theta
figure;
subplot(1,2,1)
plot(theta, P_theta)
title("Pressure in cylinder vs Angle of Crank")
xlabel("\theta (rad)");
ylabel("Pressure (Pa)");
xlim([0,max(theta)]);
xticks(0:pi/2:max(theta))
xticklabels(["0","\pi/2","\pi","3/2\pi","2\pi","5/2\pi","3\pi", ...
    "7/2\pi", "4\pi"])

%% Plot Pressure vs Volume
subplot(1,2,2)
plot(vol_theta, P_theta)
title("Pressure in Cylinder vs Volume")
xlabel("Volume (m^3)");
ylabel("Pressure (Pa)");
xlim([0.00005 0.0004]);


%% Calculate loading torque
% Calculate F from thermo equations
F_theta = P_theta*cyl_area; %Force on piston face
T_load_theta = F_theta.*AoA.*(sin(theta)-(tan(phi).*cos(theta))); % Loading torque, Nm, f(theta)


% Find the values corresponding to theta = 0 and theta = 2pi
idx_2pi = find(theta <= 2*pi, 1, 'last');

% Calculate average torque over one full crank rotation (0 to 2pi)
T_avg = (1 / (2*pi)) * trapz(theta(1:idx_2pi), T_load_theta(1:idx_2pi));
fprintf('Average Torque is %0.3f nm\n', T_avg);

% dKE = integral(T-load-T_avg, theta_0, theta_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

