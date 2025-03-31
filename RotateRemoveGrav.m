function state_corr = RotateRemoveGrav(state, curr_angle, a)
    
    measured_a_h = state(1) / cos(curr_angle);
    measured_a_v = -state(3) / cos(curr_angle) - a.g;

    state_corr = [measured_a_h measured_a_h measured_a_v];
    
end