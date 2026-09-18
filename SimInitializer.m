clc
close all

%% ---Initialize Constants--- %%
muEarth = 3.986004418e5; % Standard gravitational parameter (Earth) [km^3/s^2]

%% ---Initialize Unit Conversions--- %%
d2r = pi / 180;           % Degrees to radians
rev2r = 2*pi / (24*60^2); % Revolutions per day to radians per second