function signal_error = insertAllError(t, a)

    accel_i = a.real_accel(t);
    ang_rate_i = a.real_ang_rate(t);

    measured_a_h = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.real_ang);
    measured_a_v = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.real_ang);
    measured_g = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
    
    signal_error = [measured_a_h measured_a_h measured_a_v measured_g measured_g measured_g];

end