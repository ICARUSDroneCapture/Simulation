clear; clc; close all;

%% Set Variables
time = 0;                       % The simulation time [s]
duration = 80;                  % The simulation duration [s]
dt = 16.66/1000;                % The simulation rate [s]
simTime = 0:dt:duration;

% Function y_dot_dot = -w_n*y
omega_n = 0.1;                  % The angular velocity of the harmonic oscillator
y = [1 0 omega_n*1;          % [position, velocity, acceleration] of harmonic oscillator
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0];

state_record = zeros(5, 3);
deriv_record = zeros(5, 3);

state_int_control = zeros(length(simTime), 3);

prev_state = y(1, :);
integrator_type = 2;

%% Run the simulation
run_sim = true;
index = 1;
while run_sim == true
    if time > duration
        % break the simulation
        break;
    end

    val_dot_control = get_deriv(prev_state);

    deriv_record = insertVector(deriv_record, val_dot_control);

    state_vec_control = fdm_integrator(state_record, deriv_record, dt, integrator_type);
    
    state_record = insertVector(state_record, state_vec_control);

    state_int_control(index, :) = state_vec_control;
    prev_state = state_vec_control;
    
    % Increment sim time
    time = time + dt;
    index = index + 1;
end

y_output = state_int_control(:, 1);
y_dot_output = state_int_control(:, 2);

%% Plot the function
tiledlayout(2,1)

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


function deriv = get_deriv(state)

    y = state(2);
    y_dot = state(3);
    y_dot_dot = -0.1*y;

    deriv = [y y_dot y_dot_dot];

end

function output_vec = insertVector(originalVector, addVector)

    % Assumes the addVector is a row
    % Removes last row of originalVector
    % Addes as first row of originalVector

    output_vec = [addVector; originalVector(1:end-1, :)];

end