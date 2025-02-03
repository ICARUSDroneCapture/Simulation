%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%      Flight Dynamics Model - Harmonic Oscillator Test   %
%                                                         %
%              Principles of Flight Simulation            %
%        USED TO TEST THE STABILITY OF FDM INTEGRATORS    %
%                                                         %
%                       Jason Popich                      %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
clear
clc
tic

%% Set Variables
time = 0;                       % The simulation time [s]
duration = 80;                  % The simulation duration [s]
dt = 16.66/1000;                % The simulation rate [s]

% Function y_dot_dot = -w_n*y
omega_n = 0.1;                  % The angular velocity of the harmonic oscillator
y = [1 0 0 0 0];                % The position of the harmonic oscillator
y_dot = [0 0 0 0 0];            % The rate of change of position of the harmonic oscillator
y_dot_dot = [omega_n*y(1) 0 0 0 0];  % The rate of change of the y dot

y_output = [];
y_dot_output = [];

%% Run the simulation
run_sim = true;
index = 0;
while run_sim == true
    if time > duration
        % break the simulation
        break;
    end
    
    % Compute the rate of change
    y_dot = adjust_vector(y_dot)
    y_dot(1) = fdm_integrator(y_dot(2), y_dot_dot, dt, 3)
    
    % Compute the rate of change
    y = adjust_vector(y);
    y(1) = fdm_integrator(y(2), y_dot, dt, 3);
    
    % Update y_dot_dot
    y_dot_dot = insert_value(y_dot_dot, -omega_n*y(1))
    
    % Increment sim time
    time = time + dt;
    index = index + 1;
    
    y_output = [y_output y(1)];
    y_dot_output = [y_dot_output y_dot(1)];
end

%% Plot the function
tiledlayout(2,1)
simTime = linspace(0, time, length(y_output));

% Y Positions Plot
nexttile
plot(simTime, y_output)
hold off
title("Position Plots")
xlabel("Time [s]")
ylabel("Y")
legend('Y')

% Y Dot Positions Plot
nexttile
plot(simTime, y_dot_output)
hold off
title("Velocity Plots")
xlabel("Time [s]")
ylabel("$\dot{Y}$", 'Interpreter','latex')
legend('$\dot{Y}$', 'Interpreter','latex')


toc