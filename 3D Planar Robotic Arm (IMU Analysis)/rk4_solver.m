function [t, y] = rk4_solver(dynamics, tspan, y0, dt)
% RK4_SOLVER Solves ODEs using the fixed-step Runge-Kutta 4th-order method.
%
% Inputs:
%   dynamics - Function handle for the ODEs (dy/dt = dynamics(t, y)).
%   tspan    - Vector [t0, tf] specifying the start and end times.
%   y0       - Vector of initial conditions.
%   dt       - Fixed step size for the integration.
%
% Outputs:
%   t        - Vector of time points.
%   y        - Matrix of solution values. Each row corresponds to a time point.

    % Initialize time vector
    t = tspan(1):dt:tspan(2);
    num_steps = length(t);
    
    % Initialize solution matrix
    y = zeros(num_steps, length(y0));
    y(1, :) = y0;  % Set initial condition
    
    % RK4 Integration Loop
    for i = 1:(num_steps-1)
        ti = t(i);       % Current time
        yi = y(i, :)';   % Current solution (as a column vector)
        
        % Compute RK4 coefficients
        k1 = dt * dynamics(ti, yi);
        k2 = dt * dynamics(ti + dt/2, yi + k1/2);
        k3 = dt * dynamics(ti + dt/2, yi + k2/2);
        k4 = dt * dynamics(ti + dt, yi + k3);
        
        % Update solution
        y(i+1, :) = yi' + (k1' + 2*k2' + 2*k3' + k4') / 6;
    end
end

