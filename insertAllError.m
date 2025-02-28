function signal_error = insertAllError(t, a)

    % Uses the known deck disturbances and predefined function handles to
    accel_dot_m_v = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel(t), a.theta_err);
    accel_dot_m_h = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel(t), a.theta_err);
    gyro_dot_m = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate(t));

    signal_error = [accel_dot_m_h accel_dot_m_h accel_dot_m_v gyro_dot_m gyro_dot_m gyro_dot_m];
end