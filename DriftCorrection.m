function s_dot = DriftCorrection(a, state, signals)

    % Current states
    vel_x = state(1);
    vel_y = state(2);
    vel_z = state(3);
    theta = state(4);
    psi = state(5);
    phi = state(6);
    vel_err_accum_x = state(7);
    vel_err_accum_y = state(8);
    vel_err_accum_z = state(9);
    theta_err_accum = state(10);
    psi_err_accum = state(11);
    phi_err_accum = state(12);

    curr_state = [vel_x; vel_y; vel_z; theta; psi; phi];
    state_err_accum = [vel_err_accum_x; vel_err_accum_y; vel_err_accum_z; theta_err_accum; psi_err_accum; phi_err_accum];

    kw = a.kw; % 
    kt = a.kt; % 

    state_nominal = [0; 0; 0; 0; 0; 0];
    state_dot_nominal = [0; 0; 0; 0; 0; 0];

    accel_dot_m = signals(1:3)';
    gyro_dot_m = signals(4:6)';

    state_dot_m = [accel_dot_m; gyro_dot_m];

    state_control = curr_state - state_nominal;

    state_dot_comp = kt .* state_err_accum + kw .* state_control;

    state_dot = state_dot_m + state_dot_nominal - state_dot_comp;

    s_dot = zeros(length(state),1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end