function [t, s] = controlSystemsTesting(a)

% Simulation time
% tspan = [0 30]; % [s]

% Initial States
p0 =  a.d(a.tspan(1))+a.pr_d; % Platform position [m]
p_dot0 = a.d_dot(a.tspan(1)); % Platform velocity [m/s]
pr_err_accum0 = 0;          % Integral of error in relative position [m*s]
s0 = [p0; p_dot0; pr_err_accum0];

% Running Simulation
lastwarn('')
warnId = 'MATLAB:ode45:IntegrationTolNotMet';
op = odeset('RelTol',1e-5,'AbsTol',1e-5); % Tolerance options
% [t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);
try
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a), a.tspan,s0,op);
catch ME
    if strcmp(ME.identifier, warnId)
        % disp(beta)
        % disp(['Solver stopped due to tolerance failure at t = ', num2str(ME.stack(1).line)]);
        t = 0;
        s = 0;
    else
        % Something else has gone wrong: just re-throw the error
        throw(ME);
    end
end

end