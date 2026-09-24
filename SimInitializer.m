clc
close all

%% ---Initialize Constants--- %%
muEarth = 3.986004418e5;  % Standard gravitational parameter (Earth) [km^3/s^2]

%% ---Initialize Unit Conversions--- %%
d2r = pi / 180;           % Degrees to radians
rev2r = 2*pi / (24*60^2); % Revolutions per day to radians per second

%% ---Extract Classical Orbital Elements (COEs) from the Two Line Element (TLE)--- %%
tle = readmatrix("TLE_SharkSat.txt");

inclination = tle(2,3) * d2r;
RAAN = tle(2,4) * d2r;
eccentricity = tle(2,5);
argument_of_periapsis = tle(2,6) * d2r;
mean_anomaly = tle(2,7) * d2r;
mean_motion = tle(2,8) * rev2r;
semi_major_axis = (muEarth/mean_motion^2)^(1/3);

% Assemble COE list (7x1)
COE = [inclination; RAAN; eccentricity; argument_of_periapsis; mean_anomaly; mean_motion; semi_major_axis];

[rPQW, vPQW] = COE2RV(muEarth, COE);
orbitState0 = [rPQW; vPQW];
disp(orbitState0)
