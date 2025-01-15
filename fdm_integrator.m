function [output_value] = fdm_integrator(output, val_dot, dt, integrator_type)
%% DESCRIPTION: This function will integrate using the integrator type given
%
%% INPUT:
% output                    The vector of the previous values of the output
% val_dot                   A vector of previous dot values
% dt                        The timestep to integrate over
% integrator_type           The type of integrator to use
%                               1 - Rectangular-Euler
%                               2 - Trapezoidal Rule
%                               3 - Adams-Bashforth 2 Step Method
%                               4 - Adams-Bashforth 3 Step Method
%                               5 - Adams-Bashforth 4 Step Method
%                               6 - Adams-Bashforth 5 Step Method
%% OUTPUT: 

%% Code

% Testing single precision (equivalent to float in C++)
output = single(output);
val_dot = single(val_dot);
dt = single(dt);

switch integrator_type
    case 1
        % Rectangular-Euler
        output_value_local = output(:, 1) + dt*val_dot(:, 1);
    case 2
        % Trapezoidal with Half Step
        current_val_dot = 2.0 * val_dot(:, 1) - val_dot(:, 2);
        output_value_local = output(:, 1) + (0.5*dt*(current_val_dot + val_dot(:, 1)));
    case 3
        % Trapezoidal
        output_value_local = output(:, 1) + (0.5*dt*(val_dot(:, 1) + val_dot(:, 2)));
    case 4
        % Adams-Bashforth 2 Step
        current_val_dot = 2.0 * val_dot(:, 1) - val_dot(:, 2);
        output_value_local = output(1) + (dt*((1.5*current_val_dot) - (0.5*val_dot(:, 1))));
    case 5
        % Adams-Bashforth 3 Step
        current_val_dot = 2.0 * val_dot(:, 1) - val_dot(:, 2);
        output_value_local = output(:, 1) + ((1.0/12.0)*dt*((23.0*current_val_dot) - (16.0*val_dot(:, 1)) + (5.0*val_dot(:, 2))));
    case 6
        % Adams-Bashforth 4 Step (TODO)
        output_value_local = 0;
    case 7
        % Adams-Bashforth 5 Step (TODO)
        output_value_local = 0;
end

% Set the output as a float (single)
output_value = single(output_value_local);

end

