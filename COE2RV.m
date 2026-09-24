function [rPQW, vPQW] = COE2RV(muEarth, COE)

%% ---Extract COEs from the array--- %%
inclination = COE(1);
RAAN = COE(2);
eccentricity = COE(3);
argument_of_periapsis = COE(4);
mean_anomaly = COE(5);
mean_motion = COE(6);
semi_major_axis = COE(7);

%% ---Newton's Method--- %%

% Initialize
E = 0;
previous_E = E;
current_error = inf;
acceptable_error = 10^-6;

% Numerical Solver
while current_error > acceptable_error
    f_approx = E - eccentricity*sin(E) - mean_anomaly;
    f_deriv_approx_E = 1 - eccentricity*cos(E);

    E = previous_E - f_approx/f_deriv_approx_E;

    current_error = abs(E - previous_E);
    previous_E = E;
end
eccentric_anomaly = E;

%% ---Solve for Initial Position and Velocity--- %%

% Position
rP = semi_major_axis*(cos(E) - eccentricity);
rQ = semi_major_axis*sqrt(1-eccentricity^2) * sin(E);
rW = 0;
rPQW = [rP; rQ; rW];

% Velocity
r = norm(rPQW);
vP = -sqrt(muEarth*semi_major_axis) / r * sin(E);
vQ = sqrt(muEarth*semi_major_axis*(1-eccentricity^2)) / r * cos(E);
vW = 0;
vPQW = [vP; vQ; vW];

end