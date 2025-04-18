function signal_error = insertAllError(t, a)

    a_I = [a.real_accel_xI(t); a.real_accel_yI(t); a.real_accel_zI(t)];
    ang_rate_i = [a.theta_dot(t); a.phi_dot(t); a.psi_dot(t)];

    accel_S = a.measured_accel_3D(a, t, a_I);
    accel_S(3) = accel_S(3) + a.g;
    gyro = a.measured_gyro_3D(a, t, ang_rate_i);

    signal_error = [accel_S; gyro]';
end