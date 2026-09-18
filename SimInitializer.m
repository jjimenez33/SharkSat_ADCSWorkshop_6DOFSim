clc
close all

%% ---Initialize Constants--- %%
muEarth = 3.986004418e5; % Standard gravitational parameter (Earth) [km^3/s^2]

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

% output_data = sim('ADCSWorkshop2026_6DOFSim.slx');