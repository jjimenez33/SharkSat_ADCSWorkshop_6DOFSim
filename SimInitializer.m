clc
close all

%% ---Initialize Constants--- %%
muEarth = 3.986004418e5; % Standard gravitational parameter (Earth) [km^3/s^2]
Re      = 6378.137; % Earth's equatorial radius [km]

%% ---Initialize Unit Conversions--- %%
d2r = pi / 180;           % Degrees to radians
rev2r = 2*pi / (24*60^2); % Revolutions per day to radians per second

%% ---Extract Classical Orbital Elements (COEs) from Two Line Element (TLE) %%
tle = readmatrix("TLE_SharkSat.txt");

mean_motion = tle(2,8) * rev2r;                    % [radians/sec]
semi_major_axis = (muEarth / mean_motion^2)^(1/3); % [km]
eccentricity = tle(2,5);                           % [dimensionless]
inclination = tle(2,3)*d2r;                        % [radians]
RAAND = tle(2,4)*d2r;                              % [radians]
argument_of_periapsis = tle(2,6)*d2r;              % [radians]
mean_anomaly = tle(2,7) * d2r;                     % [radians]

% Assemble COE list Size: [7x1]
COE = [mean_motion, semi_major_axis, eccentricity, inclination, RAAND, argument_of_periapsis, mean_anomaly];

%% ---Obtain initial state vector [position; velocity] given COEs---%%
[r0, v0] = COE2RV(muEarth, COE);
orbitState0 = [r0;
               v0];

out = sim('ADCSWorkshop2026_6DOFSim.slx');

%% ---Orbit Plots--- %%

% Get time
t = out.X_Orbit_ECI.Time;

% Earth Centered Inertial %

% Get the 6xN data
X_Orbit_ECI = out.X_Orbit_ECI.Data;
% Position
rx_eci = squeeze(X_Orbit_ECI(1,1,:));
ry_eci = squeeze(X_Orbit_ECI(2,1,:));
rz_eci = squeeze(X_Orbit_ECI(3,1,:));

figure
plot(t, rx_eci)
hold on
plot(t, ry_eci)
plot(t, rz_eci)
grid on
xlabel('Time (s)')
ylabel('Position (km)')
legend('r_x', 'r_y', 'r_z')
title('ECI Position vs Time')

% Velocity
vx_eci = squeeze(X_Orbit_ECI(4,1,:));
vy_eci = squeeze(X_Orbit_ECI(5,1,:));
vz_eci = squeeze(X_Orbit_ECI(6,1,:));

figure
plot(t, vx_eci)
hold on
plot(t, vy_eci)
plot(t, vz_eci)

grid on
xlabel('Time (s)')
ylabel('Velocity (km/s)')
legend('v_x','v_y','v_z')
title('ECI Velocity vs Time')

% Earth Centered Earth Fixed %

% Get Position in ECEF
Position_ECEF = out.Position_ECEF.Data;
% Position
rx_ecef = Position_ECEF(:,1);
ry_ecef = Position_ECEF(:,2);
rz_ecef = Position_ECEF(:,3);

figure
hold on
grid on
axis equal

[Xe, Ye, Ze] = sphere(50);
globe = surf(Re*Xe, Re*Ye, -Re*Ze, 'EdgeColor', 'none');
[cdata, map] = imread('topo.png');
if ~isempty(map)
    cdata = ind2rgb(cdata, map);
end
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1)

plot3(rx_ecef, ry_ecef, rz_ecef, 'y', 'LineWidth', 2.5)
plot3(rx_ecef(1), ry_ecef(1), rz_ecef(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8)

xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)')
title('Satellite Orbit in ECEF')
view(45, 30)

