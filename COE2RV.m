function [r0, v0] = COE2RV(muEarth, COE)

%% ---Extract COEs from array--- %%
mean_motion       = COE(1);
semi_major_axis   = COE(2);
eccentricity      = COE(3);
inclination       = COE(4);
RAAN              = COE(5);
argumentPerigee   = COE(6);
mean_anomaly      = COE(7);

%% ---Newton's Method--- %%

% Initialize
E = 0;
previous_E = E;
current_error = inf;
acceptable_error = 10^-6;

% Numerical Solver
while current_error > acceptable_error
    f_approx_E       = E - eccentricity*sin(E) - mean_anomaly;
    f_deriv_approx_E = 1 - eccentricity*cos(E);

    E = previous_E - f_approx_E / f_deriv_approx_E;

    current_error = abs(E - previous_E);
    previous_E = E;
end
eccentric_anomaly = E;

%% ---Solve for Initial State Vector--- %%

% DOUBLE CHECK THESE
rP = semi_major_axis * (cos(E) - eccentricity);
rQ = semi_major_axis * sqrt(1-eccentricity^2) * sin(eccentric_anomaly);
rW = 0;
rPQW = [rP;
        rQ;
        rW];

r = semi_major_axis * (1 -eccentricity*cos(E));

vP = -sqrt(muEarth*semi_major_axis) / r * sin(E);
vQ = sqrt(muEarth*semi_major_axis*(1-eccentricity^2)) / r * cos(E);
vW = 0;
vPQW = [vP;
        vQ;
        vW];

% Rotation Matrix from Perifocal to Earth Centered Inertial Frame
R11 = cos(RAAN)*cos(argumentPerigee) - sin(RAAN)*sin(argumentPerigee)*cos(inclination);
R12 = -cos(RAAN)*sin(argumentPerigee) - sin(RAAN)*cos(argumentPerigee)*cos(inclination);
R13 = sin(RAAN)*sin(inclination);
R21 = sin(RAAN)*cos(argumentPerigee) + cos(RAAN)*sin(argumentPerigee)*cos(inclination);
R22 = -sin(RAAN)*sin(argumentPerigee) + cos(RAAN)*cos(argumentPerigee)*cos(inclination);
R23 = -cos(RAAN)*sin(inclination);
R31 = sin(argumentPerigee)*sin(inclination);
R32 = cos(argumentPerigee)*sin(inclination);
R33 = cos(inclination);

R_PQW2IJK = [R11 R12 R13;
             R21 R22 R23;
             R31 R32 R33];

% Rotate to ECI Frame
r_eci = R_PQW2IJK*rPQW;
v_eci = R_PQW2IJK*vPQW;

r0 = r_eci;
v0 = v_eci;


end





