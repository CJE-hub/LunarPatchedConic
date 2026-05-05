%% Lunar Patched Conic Trajectory (w/ Plots)

clear all; clc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Conversions

deg    = pi/ 180;    %...Deg --> Rad
meters = 1000;       %...Km  --> M
hr     = 3600;       %...Sec --> Hr
g0 = 9.80665 / 1000; %...Standard gravity (km/s^2)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Constants

D          = 3.844 * 10^5;        %...Mean Earth / Moon Distance
R_Moon_SOI = 6.63 * 10^4;         %...Radius of Lunar SOI
r_moon     = 1737.370343;         %...Radius of Moon
r_earth    = 6378.1363;           %...Radius of Earth
omega_moon = 2.66 * 10^(-6);      %...Lunar Angular Velocity
V_m        = 1.018;               %...Magnitude of Geocendric Lunar Velocity
mu_earth   = 3.986004415 * 10^5;  %...Earth Gravational Constant 
mu_moon    = 4.902843788 * 10^3 ; %...Lunar Gravational Constant

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Input Geocentric Initial Conditions (Assuming Circular)

zp  = 318.9073;      %...Initial Orbit Periapsis Altitude
za  = 35000;         %...Initial Orbit Apoapsis Altitude
inc_earth  = 28.5;   %...Earth Parking Orbit Inclination
raan_earth = 0;      %...Earth Parking Orbit RAAN
w_earth    = 0;      %...Earth Parking Orbit Arg of Periapsis
r_0 = zp + r_earth;  %...Magnitude of Position Vector at Departure
v_0 = 10.84618;      %...Magnitude of Velocity Vector at Departure
FPA = 0;             %...Flight Path Angle

%...Target Orbit Parameter (Moon)
target_perilune = 100;  %...Target Altitude of Perilune (km)
target_apolune  = 100;  %...Target Altitude of Apolune (km)
inc_moon  = 0;          %...Target Lunar Inclination (Deg)
raan_moon = 0;          %...Target Lunar RAAN (Deg)

%...Spacecraft Parameters
Isp = 300;               %...Specific Impulse (sec)
m_initial = 5000;        %...Inital S/C Mass in LEO (kg)

%...Total Allowable Mission Time 
max_mission_time = 5;  %..Days

%...End Input
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Orbital Parameters of Initial Earth Parking Orbit

rp_earth = r_0;
ra_earth = za + r_earth;
ra_leo = r_earth + za;
a_leo  = (rp_earth + ra_earth) / 2;
e_leo  = (ra_earth - rp_earth) / (ra_earth + rp_earth);
Vp_leo = sqrt(mu_earth * ((2 / rp_earth) - (1 / a_leo)));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Orbital Parameters of Geocentric Lunar Transfer Ellipse

E   = ((v_0^2) / 2)  - (mu_earth / r_0); %...Specific Mechanical Energy
h   = r_0*v_0*cos(FPA);                  %...Specific Angular Momentum
p   = (h^2) / mu_earth;                  %...Semi-Parameter
a   = - mu_earth / (2 * E);              %...Semi-major Axis
e   = sqrt(1 - (p / a));                 %...Eccentricity

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Arrival Angle Optimization Search 

target_rp = target_perilune + r_moon;    % Target Perilune Radius (km)

%...Define the search space for Arrival Angle (Degrees)
AA_search = linspace(10, 80, 5000);      %...10-80 degrees with 5000 steps
rp_search = zeros(size(AA_search));

for i = 1:length(AA_search)
    AA_test = AA_search(i);
    
    %...Test Geocentric SOI Arrival
    r_1_test   = sqrt((D^2) + (R_Moon_SOI^2) ...
        - (2 * D * R_Moon_SOI * cos(AA_test * deg)));
    V_1_test   = sqrt(2 * ( E + (mu_earth / r_1_test)));
    FPA_a_test = acos(h / (r_1_test * V_1_test));
    phi_1_test = asin((R_Moon_SOI / r_1_test)*sin(AA_test * deg));
    
    %...Test Selenocentric Conversion
    V_2_test = sqrt((V_1_test^2) + (V_m^2) - ...
        (2 * V_1_test * V_m * cos((FPA_a_test - phi_1_test))));
    e_2_test = asin(((V_m / V_2_test) * cos(AA_test*deg)) - ...
        ((V_1_test / V_2_test) ...
        * cos(((AA_test*deg) + phi_1_test - FPA_a_test))));
    
    %...Test Lunar Orbital Parameters
    E_moon_test = ((V_2_test^2) / 2)  - (mu_moon / R_Moon_SOI);
    h_moon_test = R_Moon_SOI * V_2_test * sin(e_2_test);
    p_moon_test = (h_moon_test^2) / mu_moon;
    a_moon_test = - mu_moon / (2 * E_moon_test);
    e_moon_test = sqrt(1 - (p_moon_test / a_moon_test));
    
    %...Store Successful Perilune Radius
    rp_search(i) = p_moon_test / (1 + e_moon_test);
end

%...Search Closest Result
[min_miss, best_idx] = min(abs(rp_search - target_rp));

%...Redefine with Optimal Result
AA = AA_search(best_idx); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Geocentric Conditions at Lunar Sphere of Influence

%...Magnitude of Position Vector @ Arrival
r_1   = sqrt( (D^2) + (R_Moon_SOI^2) ...
    - (2 * D * R_Moon_SOI * cos(AA * deg)));
%...Magnitude of Velocity Vector @ Arrival
V_1   = sqrt(2 * ( E + (mu_earth / r_1)));
%...Flight Path Anlge @ Arrival
FPA_a = acos(h / (r_1 * V_1));
%...Phase Angle @ Arrival
phi_1 = asin((R_Moon_SOI / r_1)*sin(AA * deg));

%...DV for Trans-Lunar Injection
DV_translunar =  v_0 - Vp_leo;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Time of Flight to Lunar Sphere of Influence 

%...Clamp arguments between -1 and 1 to prevent floating-point complex errors
arg0 = max(min((p - r_0) / (r_0 * e), 1), -1);
arg1 = max(min((p - r_1) / (r_1 * e), 1), -1);

%...True Anomaly @ Departure & Arrival
TA_0 = acos(arg0);
TA_1 = acos(arg1);

%...True Anomaly @ Departure
TA_0 = acos((p - r_0) / (r_0 * e));
%...True Anomaly @ Arrival
TA_1 =acos((p - r_1) / (r_1 * e));
%...Eccentric Anomaly @ Departure
E_0 = 2 * atan2(sqrt(1-e)*sin(TA_0/2), sqrt(1+e)*cos(TA_0/2));
%...Eccentric Anomaly @ Arrival
E_1 = 2 * atan2(sqrt(1-e)*sin(TA_1/2), sqrt(1+e)*cos(TA_1/2));
%...Time of Flight
dt = sqrt((a^3) / mu_earth) * ( (E_1 - (e*sin(E_1)) ...
    -  (E_0 - (e*sin(E_0)))));
%...Phase Angle @ Departure
phi_0 = TA_1 - TA_0 - phi_1 - (omega_moon*dt);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Convert Geocentric LSI Arrival Conditions to Selenocentric Conditions

R_2 = R_Moon_SOI;                               %...Magnitude of Position Vector
V_2 = sqrt((V_1^2) + (V_m^2) ...
    - (2 * V_1 * V_m * cos((FPA_a - phi_1 )))); %...Magnitude of Velocity Vector    
e_2 = asin( ((V_m / V_2)*cos(AA*deg))  ...      %...Zenith Angle to Lunar Center
    -   ( (V_1 / V_2)*cos(((AA*deg) + phi_1 - FPA_a))));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Lunar Orbital Elements

E_moon = ((V_2^2) / 2)  - (mu_moon / R_2);  %...Specific Mechanical Energy
h_moon = R_2*V_2*sin(e_2);                  %...Specific Angular Momentum
p_moon = (h_moon^2) / mu_moon;              %...Semi-Parameter
a_moon = - mu_moon / (2 * E_moon);          %...Semi-major Axis
e_moon = sqrt(1 - (p_moon / a_moon));       %...Eccentricity

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Lunar Orbit Insertion Parameters (Elliptical Target)

rp_moon = p_moon / (1 + e_moon);                   %...Radius of Perilune (Arrival)
zp_moon = rp_moon - r_moon;                        %...Altitude of Perilune
Vp_moon = sqrt((V_2^2) + (2 * (mu_moon/rp_moon))); %...Velocity at Perilune (Arrival)

%...Target Elliptical Orbit
ra_moon_target = target_apolune + r_moon;
a_moon_target = (rp_moon + ra_moon_target) / 2;
e_moon_target = (ra_moon_target - rp_moon) / (ra_moon_target + rp_moon);

%...Required Velocity at Perilune for Target Apoapsis
Vp_target = sqrt(mu_moon * ((2 / rp_moon) - (1 / a_moon_target)));

%...DV Lunar Burn
DV_lunar = Vp_moon - Vp_target;                       

%...Total DV
DV_tot = abs(DV_translunar) + abs(DV_lunar);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Total Mission Time

%...Time Vector for One Full Orbit Period
Period_moon = 2 * pi * sqrt(rp_moon^3 / mu_moon);

%...Calculate Total Mission Time (1 LEO Phasing + Transfer + 1 LLO Phasing)
T_LEO = 2 * pi * sqrt(a_leo^3 / mu_earth);
total_time_sec = T_LEO + dt + Period_moon;
total_time_days = total_time_sec / 86400;

%...ComputeTime Constraint
if total_time_days <= max_mission_time
    time_status = 'GO';
else
    time_status = 'NO GO';
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Trajectory Data for Plotting (Keplers Equation) 

theta_circle = linspace(0, 2*pi, 200);

%...Initial Earth Parking Orbit Array
Period_leo = 2 * pi * sqrt(a_leo^3 / mu_earth);
t_leo = linspace(0, Period_leo, 200);
n_leo = sqrt(mu_earth / a_leo^3);

x_leo_pqw = zeros(size(t_leo));
y_leo_pqw = zeros(size(t_leo));

%...Propagate the true anomaly over one full period
for i = 1:length(t_leo)
    M_i = n_leo * t_leo(i);    
    E_i = kepler_E(e_leo, M_i);    
    nu_i = 2 * atan2(sqrt(1+e_leo)*sin(E_i/2), sqrt(1-e_leo)*cos(E_i/2));
    r_i = a_leo * (1 - e_leo * cos(E_i));     
    x_leo_pqw(i) = r_i * cos(nu_i);
    y_leo_pqw(i) = r_i * sin(nu_i);
end

%...Rotate to 3D Geocentric Frame(Euler Rotation)
inc_e_rad  = inc_earth * deg;
raan_e_rad = raan_earth * deg;
w_e_rad    = w_earth * deg;

x_leo = (cos(raan_e_rad)*cos(w_e_rad) - sin(raan_e_rad)*sin(w_e_rad)*cos(inc_e_rad)) .* x_leo_pqw ...
      + (-cos(raan_e_rad)*sin(w_e_rad) - sin(raan_e_rad)*cos(w_e_rad)*cos(inc_e_rad)).* y_leo_pqw;
y_leo = (sin(raan_e_rad)*cos(w_e_rad) + cos(raan_e_rad)*sin(w_e_rad)*cos(inc_e_rad)) .* x_leo_pqw ...
      + (-sin(raan_e_rad)*sin(w_e_rad) + cos(raan_e_rad)*cos(w_e_rad)*cos(inc_e_rad)) .* y_leo_pqw;
z_leo = (sin(w_e_rad)*sin(inc_e_rad)) .* x_leo_pqw ...
      + (cos(w_e_rad)*sin(inc_e_rad)) .* y_leo_pqw;

%...Geocentric Transfer Trajectory Array 
t_vec = linspace(0, dt, 500);          % Time of Transfer
n     = sqrt(mu_earth / a^3);          % Mean Motion of the Transfer Orbit 
M_0   = E_0 - e * sin(E_0);            % Mean Anomaly @ Departure

x_traj_pqw = zeros(size(t_vec));          
y_traj_pqw = zeros(size(t_vec));

%...2D Perifocal Plane 
for i = 1:length(t_vec)
    M_i = M_0 + n * t_vec(i);    
    E_i = kepler_E(e, M_i);    
    nu_i = 2 * atan2(sqrt(1+e)*sin(E_i/2), sqrt(1-e)*cos(E_i/2));
    r_i = a * (1 - e * cos(E_i));     
    x_traj_pqw(i) = r_i * cos(nu_i);
    y_traj_pqw(i) = r_i * sin(nu_i);
end

%...Rotate Geocentric Transfer to match 3D Earth Parking Orbit
x_traj = (cos(raan_e_rad)*cos(w_e_rad) - sin(raan_e_rad)*sin(w_e_rad)*cos(inc_e_rad)) .* x_traj_pqw ...
       + (-cos(raan_e_rad)*sin(w_e_rad) - sin(raan_e_rad)*cos(w_e_rad)*cos(inc_e_rad)) .* y_traj_pqw;
y_traj = (sin(raan_e_rad)*cos(w_e_rad) + cos(raan_e_rad)*sin(w_e_rad)*cos(inc_e_rad)) .* x_traj_pqw ...
       + (-sin(raan_e_rad)*sin(w_e_rad) + cos(raan_e_rad)*cos(w_e_rad)*cos(inc_e_rad)) .* y_traj_pqw;
z_traj = (sin(w_e_rad)*sin(inc_e_rad)) .* x_traj_pqw ...
       + (cos(w_e_rad)*sin(inc_e_rad)) .* y_traj_pqw;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Moon's Inertial Syncing & Phase Angle Calculation

%...Compute Exact Angle where S/C Hits Lunar SOI
theta_arrival = atan2(y_traj(end), x_traj(end));

%...Moon MUST be at Exact Angle (t = dt)
x_moon1 = D * cos(theta_arrival);
y_moon1 = D * sin(theta_arrival);

%...Calculate where the Moon was at t = 0 (Launch)
theta_moon_launch = theta_arrival - (omega_moon * dt);


%...(0-360 degrees) for dashboard table
phi_0_deg = mod(theta_moon_launch * 180/pi, 360); 

%...Moon's Full Orbit Array 
theta_orb = linspace(0, 2*pi, 300);
x_orb = D * cos(theta_orb);
y_orb = D * sin(theta_orb);
z_orb = zeros(size(x_orb));

%...Lunar SOI Boundary at Arrival Plot
x_soi = x_moon1 + R_Moon_SOI * cos(theta_circle);
y_soi = y_moon1 + R_Moon_SOI * sin(theta_circle);

%...Selenocentric (Moon-Centered) Array
x_m = r_moon * cos(theta_circle);
y_m = r_moon * sin(theta_circle);
x_soi_moon = R_Moon_SOI * cos(theta_circle);
y_soi_moon = R_Moon_SOI * sin(theta_circle);

%...Convert Target Angles to Radians
inc_rad = inc_moon * deg;
raan_rad = raan_moon * deg;

%...Hyperbolic Arrival Trajectory Array (Rotated into 3D Target Plane)
theta_SOI = acos((p_moon / R_Moon_SOI - 1) / e_moon);
theta_hyp = linspace(-theta_SOI, 0, 300);
r_hyp = p_moon ./ (1 + e_moon * cos(theta_hyp));

%...2D Perifocal Plane
x_hyp_pqw = r_hyp .* cos(theta_hyp);
y_hyp_pqw = r_hyp .* sin(theta_hyp);

%...Rotate to 3D Selenocentric Frame (Euler Rotation)
x_hyp = (cos(raan_rad) .* x_hyp_pqw) - (sin(raan_rad) .* cos(inc_rad) .* y_hyp_pqw);
y_hyp = (sin(raan_rad) .* x_hyp_pqw) + (cos(raan_rad) .* cos(inc_rad) .* y_hyp_pqw);
z_hyp = (sin(inc_rad) .* y_hyp_pqw);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Final Lunar Parking Orbit Array (Rotated into 3D Target Plane)

%...Define COEs for Parking Orbit at the exact moment of LOI
h_park = rp_moon * Vp_target;    % Angular momentum of orbit
e_park = e_moon_target;          % orbit eccentricity
w_park = 0;                      % Argument of periapsis (0 for circular)
TA_loi = 0;                      % True Anomaly is 0 at the LOI burn point

coe_park = [h_park, e_park, raan_rad, inc_rad, w_park, TA_loi];

%...Extract Cartesian State Vector [R, V] 
[r_loi, v_loi] = sv_from_coe(coe_park, mu_moon);

%...Single 6x1 for ODE45
Y0 = [r_loi, v_loi]'; 

%...Numerical Integration (ODE45 Propagator)

%...Time vector for one full elliptical orbit period
Period_moon_target = 2 * pi * sqrt(a_moon_target^3 / mu_moon);
t_park = linspace(0, Period_moon_target, 2000); % Bumped to 2000 for a smoother curve

%...Two-Body equations of motion as an inline anonymous function
two_body_ode = @(t, Y) [Y(4:6); (-mu_moon / norm(Y(1:3))^3) * Y(1:3)];

%...Propagate
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, Y_out] = ode45(two_body_ode, t_park, Y0, options);

%...Extract the propagated trajectory data (Plotting)
x_park = Y_out(:, 1)';
y_park = Y_out(:, 2)';
z_park = Y_out(:, 3)';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Mass Fractions and Propellant Required

%...TLI Burn
dV_tli = abs(DV_translunar);
m_after_tli = m_initial * exp(-dV_tli / (Isp * g0));
m_prop_tli = m_initial - m_after_tli;

%...LOI Burn
dV_loi = abs(DV_lunar);
m_after_loi = m_after_tli * exp(-dV_loi / (Isp * g0));
m_prop_loi = m_after_tli - m_after_loi;

%...Total Propellant
m_prop_total = m_prop_tli + m_prop_loi;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%...Plots, Figures and Data

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Create UI Panel

figure('Position', [100, 100, 1400, 900], 'Name', 'Lunar Mission Data Screen' ...
    , 'NumberTitle', 'off');
hold on;

%...Define uipanels for left (plots) and right (data) sections
left = uipanel('Parent', gcf, 'Position', [0.01 0.01 0.60 0.98], 'Title', ...
    'Mission Trajectory Plots', 'BackgroundColor', [1 1 1]);
right = uipanel('Parent', gcf, 'Position', [0.62 0.01 0.37 0.98], 'Title', ...
    'Mission Data & ΔV Budget', 'BackgroundColor', [1 1 1]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...LEFT SECTION: Subplots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%...Helper Sphere Data
[X_sphere, Y_sphere, Z_sphere] = sphere(30);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...SUBPLOT 1: Earth Departure & TLI
hax1a = axes('Parent', left, 'Position', [0.02, 0.55, 0.44, 0.40]); 
hold(hax1a, 'on'); grid(hax1a, 'on'); axis(hax1a, 'equal'); 
view(hax1a, [ 50.6696, 23.9617]);
title(hax1a, 'Earth Departure & TLI Burn');
xlabel(hax1a, 'X (km)'); ylabel(hax1a, 'Y (km)'); zlabel(hax1a, 'Z (km)');

%...Draw the Geocentric Equatorial Plane (Z=0)
plane_radius = 40000; % Extended past the 35,000km apoapsis for contrast
theta_plane = linspace(0, 2*pi, 100);
fill3(hax1a, plane_radius * cos(theta_plane), plane_radius * sin(theta_plane), zeros(1, 100), ...
    [0.8 0.8 0.8], 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineStyle', '--', 'DisplayName', 'Equatorial Plane');


%...Draw the Earth and Trajectories
surf(hax1a, X_sphere * r_earth, Y_sphere * r_earth, Z_sphere * r_earth, ...
    'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'Earth');
plot3(hax1a, x_leo, y_leo, z_leo, 'c-', 'LineWidth', 1.5, ...
    'DisplayName', 'Initial Elliptical Orbit');
plot3(hax1a, x_traj, y_traj, z_traj, 'b-', 'LineWidth', 2, 'DisplayName', ...
    'Outgoing Trajectory');

%...Dynamically calculate the 3D vector (TLI Arrow)
tli_x = x_traj(1); tli_y = y_traj(1); tli_z = z_traj(1);
dir_x_tli = x_traj(2) - x_traj(1);
dir_y_tli = y_traj(2) - y_traj(1);
dir_z_tli = z_traj(2) - z_traj(1);
mag_tli = norm([dir_x_tli, dir_y_tli, dir_z_tli]);

plot3(hax1a, tli_x, tli_y, tli_z, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'TLI Point');
quiver3(hax1a, tli_x, tli_y, tli_z, (dir_x_tli/mag_tli)*15000, (dir_y_tli/mag_tli)*15000, (dir_z_tli/mag_tli)*15000, ...
    'r', 'LineWidth', 2, 'MaxHeadSize', 2, 'DisplayName', 'TLI \DeltaV');

axis(hax1a, [-45000 45000 -45000 45000 -25000 25000]);
camlight(hax1a, 'headlight'); lighting(hax1a, 'gouraud'); legend(hax1a, 'Location', 'best');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...SUBPLOT 2: Full Geocentric Transfer
hax1b = axes('Parent', left, 'Position', [0.52, 0.55, 0.44, 0.40]); 
hold(hax1b, 'on'); grid(hax1b, 'on'); axis(hax1b, 'equal'); 
view(hax1b, [-90, 90]);
title(hax1b, 'Earth-Moon Transfer (To Scale)');
xlabel(hax1b, 'X (km)'); ylabel(hax1b, 'Y (km)'); zlabel(hax1b, 'Z (km)');

surf(hax1b, X_sphere * r_earth, Y_sphere * r_earth, Z_sphere * r_earth, ...
    'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'Earth');
plot3(hax1b, x_traj, y_traj, z_traj, 'b-', 'LineWidth', 2, 'DisplayName', 'Transfer Trajectory');
plot3(hax1b, x_orb, y_orb, z_orb, 'k--', 'Color', [0 0 0 0.5], 'DisplayName', 'Moons Orbit');

%...Use the Synced Moon Coordinates
surf(hax1b, (X_sphere * r_moon) + x_moon1, (Y_sphere * r_moon) + y_moon1, Z_sphere * r_moon, ...
    'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'Moon');
plot3(hax1b, x_soi, y_soi, zeros(size(x_soi)), 'y:', 'LineWidth', 2, 'DisplayName', 'Lunar SOI');

camlight(hax1b, 'headlight'); lighting(hax1b, 'gouraud'); 
legend(hax1b, 'Location', 'best');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...SUBPLOT 3: Lunar Arrival 2D

hax2 = axes('Parent', left, 'Position', [0.02, 0.05, 0.44, 0.40]); 
hold(hax2, 'on'); grid(hax2, 'on'); axis(hax2, 'equal');
title(hax2, 'Moon-Centered Arrival');
xlabel(hax2, 'X (km)'); ylabel(hax2, 'Y (km)');

fill(hax2, x_m, y_m, [0.7 0.7 0.7], 'DisplayName', 'Moon');
plot(hax2, x_soi_moon, y_soi_moon, 'y:', 'LineWidth', 1.5, 'DisplayName', 'Lunar SOI');
plot(hax2, x_hyp, y_hyp, 'cyan', 'LineWidth', 2, 'DisplayName', 'Hyperbolic Trajectory');
plot(hax2, x_park, y_park, 'b-', 'LineWidth', 2, 'DisplayName', 'Final Lunar Orbit');
plot(hax2, x_hyp(end), y_hyp(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'LOI Burn');

%...Scale to Fit SOI
bound_2d = max(1.1 * R_Moon_SOI, 1.1 * target_apolune);
axis(hax2, [-bound_2d, bound_2d, -bound_2d, bound_2d]);
%axis(hax2, [-1.1*R_Moon_SOI, 1.1*R_Moon_SOI, -1.1*R_Moon_SOI,
%1.1*R_Moon_SOI]); %...For Lower Alt Orbit
% Zoom in closely to show the curve of the hyperbola and the LLO
%axis(hax2, [-10000, 10000, -10000, 10000]);
legend(hax2, 'Location', 'best');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...SUBPLOT 4: Lunar Arrival & LOI Burn

hax3 = axes('Parent', left, 'Position', [0.52, 0.05, 0.44, 0.40]); 
hold(hax3, 'on'); grid(hax3, 'on'); axis(hax3, 'equal'); 
view(hax3, [ 124.4411, 36.2517]); % Matches the Earth departure viewing angle
title(hax3, 'Moon Arrival & LOI Burn');
xlabel(hax3, 'X (km)'); ylabel(hax3, 'Y (km)'); zlabel(hax3, 'Z (km)');

surf(hax3, X_sphere * r_moon, Y_sphere * r_moon, Z_sphere * r_moon, ...
    'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'DisplayName', 'Moon');
%plot3(hax3, x_park, y_park, z_park, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Final LLO');
%plot3(hax3, x_hyp, y_hyp, z_hyp, 'cyan', 'LineWidth', 2, 'DisplayName', 'Incoming Hyperbola');
plot3(hax3, x_park, y_park, z_park, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Final Lunar Orbit');
plot3(hax3, x_hyp, y_hyp, z_hyp, 'cyan', 'LineWidth', 2, 'DisplayName', 'Incoming Hyperbola');
%...Find 3D coordinates of the LOI point
loi_x = x_hyp(end); 
loi_y = y_hyp(end); 
loi_z = z_hyp(end);
plot3(hax3, loi_x, loi_y, loi_z, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'LOI Point');

dir_x = x_hyp(end-1) - x_hyp(end);
dir_y = y_hyp(end-1) - y_hyp(end);
dir_z = z_hyp(end-1) - z_hyp(end);

%...Normalize the vector and scale
vec_mag = sqrt(dir_x^2 + dir_y^2 + dir_z^2);
%quiver3(hax3, loi_x, loi_y, loi_z, (dir_x/vec_mag)*2500, (dir_y/vec_mag)*2500, (dir_z/vec_mag)*2500, ...
%    'r', 'LineWidth', 2, 'MaxHeadSize', 2, 'DisplayName', 'LOI \DeltaV');
%axis(hax3, [-6000 6000 -6000 6000 -5000 5000]); % Zoomed in closely to the Moon
bound_3d = max(6000, 1.1 * target_apolune);
axis(hax3, [-bound_3d bound_3d -bound_3d bound_3d -bound_3d bound_3d]);
camlight(hax3, 'headlight'); lighting(hax3, 'gouraud'); legend(hax3, 'Location', 'best');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...RIGHT SECTION: Tables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Table 1: Initial Conditions & Time of Flight

Table1_ColNames = {'Parameter', 'Value', 'Units'};
Table1_Data = {...
    'Initial Parking Periapsis Alt', sprintf('%.4f', zp), 'km'; ...
    'Required Departure Velocity', sprintf('%.5f', v_0), 'km/sec'; ...
    'Initial Flight Path Angle', sprintf('%.2f', FPA * 180/pi), 'deg'; ...
    'Required Lunar Phase Angle at Launch', sprintf('%.2f', phi_0_deg), 'deg'; ...
    'Earth Phasing Orbit Period (1 LEO)', sprintf('%.2f', T_LEO / 3600), 'Hours'; ...
    'Earth-Moon Transfer Time of Flight', sprintf('%.2f', dt / 3600), 'Hours'; ...
    'Lunar Phasing Orbit Period (1 LLO)', sprintf('%.2f', Period_moon / 3600), 'Hours'; ...
    'Total Mission Duration', sprintf('%.2f', total_time_days), 'Days'; ...
    'Max Mission Time Constraint', sprintf('%.1f', max_mission_time), 'Days'; ...
    'Mission Time Status', time_status, '---' ...
};

%...Shaping UI
t1_parent = uipanel('Parent', right, 'Position', [0.02, 0.69, 0.96, 0.26], 'Title', ...
    'Table 1: Initial Conditions & Constants');
uitable('Parent', t1_parent, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', Table1_Data, 'ColumnName', Table1_ColNames, ... 
    'ColumnWidth', {220, 100, 80}, 'RowName', []);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Table 2: Optimization Parameters

Table2_ColNames = {'Optimization Metric', 'Value', 'Units'};
Table2_Data = {...
    'Target Perilune Altitude', sprintf('%.2f', target_perilune), 'km'; ...
    'Target Perilune Radius', sprintf('%.2f', target_perilune + r_moon), 'km'; ...
    'Optimal Arrival Angle (Calculated)', sprintf('%.4f', AA), 'deg'; ...
    'Achieved Perilune Altitude', sprintf('%.4f', zp_moon), 'km' ...
};

%...Shaping UI
t2_parent = uipanel('Parent', right, 'Position', [0.02, 0.55, 0.96, 0.13], 'Title', ...
    'Table 2: Arrival Angle Optimization Results');
uitable('Parent', t2_parent, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', Table2_Data, 'ColumnName', Table2_ColNames, ... 
    'ColumnWidth', {220, 100, 80}, 'RowName', []);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Table 3: Lunar Trajectory Velocities

Table3_ColNames = {'#', 'Description', 'Velocity (km/sec)'};
Table3_Data = {...
    '1', 'Initial Elliptical Periapsis', sprintf('%.6f', Vp_leo); ...
    '2', 'Translunar Injection Required Velocity', sprintf('%.5f', v_0); ...
    '3', 'Geocentric Velocity at Lunar Sphere of Influence', sprintf('%.6f', V_1); ...
    '4', 'Selenocentric Velocity at Lunar SOI', sprintf('%.6f', V_2); ...
    '5', 'Velocity at Perilune', sprintf('%.6f', Vp_moon); ...
    '6', 'Target Orbit Periapsis Velocity', sprintf('%.6f', Vp_target); ...
};

%...Shaping UI
t3_parent = uipanel('Parent', right, 'Position', [0.02, 0.36, 0.96, 0.18], 'Title', ...
    'Table 3: Lunar Trajectory Velocities');
uitable('Parent', t3_parent, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', Table3_Data, 'ColumnName', Table3_ColNames, ... 
    'ColumnWidth', {25, 275, 120}, 'RowName', []);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Table 4: Sample Lunar Trajectory ΔV Budget

Table4_ColNames = {'#', 'Description', 'ΔV (km/sec)'};
Table4_Data = {...
    '1', 'Launch to Phasing Orbit', sprintf('%.6f km/sec', Vp_leo); ...
    '2', 'Translunar Injection Burn (TLI)', sprintf('%.6f km/sec', abs(DV_translunar)); ...
    '3', 'Lunar Orbit Insertion burn (LOI)', sprintf('%.6f km/sec', abs(DV_lunar)); ...
    '', 'Total ΔV Budget ', sprintf('%.6f km/sec', DV_tot) ...
};

%...Shaping UI
t4_parent = uipanel('Parent', right, 'Position', [0.02, 0.22, 0.96, 0.13], 'Title', ...
    'Table 4: Lunar Trajectory ΔV Budget');
uitable('Parent', t4_parent, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', Table4_Data, 'ColumnName', Table4_ColNames, ... 
    'ColumnWidth', {25, 275, 120}, 'RowName', []);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%...Table 5: Spacecraft Mass & Propellant Budget

Table5_ColNames = {'Mass Component', 'Value', 'Units'};
Table5_Data = {...
    'Assumed Engine Specific Impulse (Isp)', sprintf('%d', Isp), 'sec'; ...
    'Initial Mass in LEO (IMLEO)', sprintf('%.1f', m_initial), 'kg'; ...
    'Propellant Consumed for TLI', sprintf('%.1f', m_prop_tli), 'kg'; ...
    'Propellant Consumed for LOI', sprintf('%.1f', m_prop_loi), 'kg'; ...
    'Final Spacecraft Dry Mass in Lunar Orbit', sprintf('%.1f', m_after_loi), 'kg'; ...
    'Total Propellant Required', sprintf('%.1f', m_prop_total),'kg' ...
};
%...Shaping UI

t5_parent = uipanel('Parent', right, 'Position', [0.02, 0.02, 0.96, 0.19], 'Title', ...
    'Table 5: Propellant Mass Budget (Rocket Equation)');
uitable('Parent', t5_parent, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', Table5_Data, 'ColumnName', Table5_ColNames, ... 
    'ColumnWidth', {220, 100, 80}, 'RowName', []);

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%... DATA EXPORT UTILITY
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Define the name for this specific run (Change this before running!)
Run_Name = 'HEO_to_Equatorial'; 
%Run_Name = 'HEO_to_NRHO'; 

%...Package the critical mission parameters into a structure
ExportData = struct();
ExportData.Name = Run_Name;
ExportData.Isp = Isp;
ExportData.DV_TLI = abs(DV_translunar);
ExportData.DV_LOI = abs(DV_lunar);
ExportData.DV_Total = DV_tot;
ExportData.TOF_Days = total_time_days;

%...Save the structure as a .mat file in the current directory
filename = sprintf('Export_%s.mat', Run_Name);
save(filename, '-struct', 'ExportData');

fprintf('\n>>> SUCCESS: Mission data saved to %s <<<\n', filename);

%% Created by Cameron Edwards, Florida Institute of Technology (April - May 2026)