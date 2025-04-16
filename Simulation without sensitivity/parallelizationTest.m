close all; clc; clear;

SimulationParameters;

% solving the system
dt = 0.001; %[d]
time_interval = [0 80]; %seconds

initial_conditions = [0; 0; 0; 
                       0; 0; 0]; %[q1; Dq1; int_q1_err; 
                                  % pm_ddot];

% Computation time: 12.8444 hrs
acceleration_gains = 0:0.5:5;
rel_prop_gains = 0:0.2:1;
rel_deriv_gains = 0:1:10;
rel_int_gains = 0:0.1:1;
% acceleration_gains = 1;
% rel_prop_gains = 0.1;
% rel_deriv_gains = 5;
% rel_int_gains = 0.5;

totalIter = length(rel_int_gains) * length(rel_deriv_gains) * ...
            length(rel_prop_gains) * length(acceleration_gains);

% Use cell arrays for parfor compatibility
max_isolation_temp = cell(totalIter, 1);

% If single simulation run, plot response over time
plotSim = false;
if totalIter == 1
    plotSim = true;
    % Initialize variables to hold plotting data
    t_wave_temp = cell(1);
    s_wave_temp = cell(1);
    a_wave_temp = cell(1);
    isolation_wave_temp = cell(1);
end

if isempty(gcp('nocreate'))
    if plotSim
        parpool('local',1);
    else
        parpool('local');
    end
end

tic;

parfor idx_all = 1:totalIter

    [a_idx, p_idx, d_idx, i_idx] = ind2sub([length(acceleration_gains), ...
                                                          length(rel_prop_gains), ...
                                                          length(rel_deriv_gains), ...
                                                          length(rel_int_gains)], ...
                                                          idx_all);

    atemp = a;
    atemp.ka = acceleration_gains(a_idx);
    atemp.kp_c = rel_prop_gains(p_idx);
    atemp.kd_c = rel_deriv_gains(d_idx);
    atemp.ki_c = rel_int_gains(i_idx);

    MFun = @(t, y)EOM_V3(t, y, atemp);
    [t, s]= rk4_solver(MFun, time_interval, initial_conditions, dt);
    s = s';

    if plotSim
        t_wave_temp{idx_all} = t
        s_wave_temp{idx_all} = s
        a_wave_temp{idx_all} = atemp;
    end

    theta2_eval = atemp.thetad(t);
    
    zEE = 0 - atemp.l1*sin(s(1,:)+theta2_eval);
    zNotIso = 0 - atemp.l1*sin(atemp.q1_ref+theta2_eval);
    
    t_test = 20;
    isolation = calculateAverageIsolation(zEE(t > t_test), ...
                                          zNotIso(t > t_test));

    if isolation < 0
        isolation = NaN;
    end

    if plotSim
        isolation_wave_temp{idx_all} = isolation;
    end

    if ~isnan(isolation)
        max_isolation_temp{idx_all} = isolation;
    else
        max_isolation_temp{idx_all} = NaN;
    end
end

if plotSim
    t_wave = t_wave_temp{1};
    s_wave = s_wave_temp{1};
    a_wave = a_wave_temp{1};
    isolation_wave = isolation_wave_temp{1};
    max_isolation = max_isolation_temp{1};
    % plotResponse(t_wave, s_wave, betas, isolation_wave, a_wave, max_isolation)
end

%% Reconstruct max_isolation after parfor
max_isolation = NaN(length(acceleration_gains), length(rel_prop_gains), ...
                    length(rel_deriv_gains), length(rel_int_gains));
for idx_all = 1:totalIter
    [a_idx, p_idx, d_idx, i_idx, rg_idx, B_idx, n_idx] = ind2sub([length(acceleration_gains), ...
                                                          length(rel_prop_gains), ...
                                                          length(rel_deriv_gains), ...
                                                          length(rel_int_gains)], ...
                                                          idx_all);
    if max_isolation_temp{idx_all}
    max_isolation(a_idx, p_idx, d_idx, i_idx) = max_isolation_temp{idx_all};
    end
end

toc

% Map indices to actual axis values
x1 = acceleration_gains;
y1 = rel_prop_gains;
z1 = rel_deriv_gains;
a1 = rel_int_gains;

[M,I] = min(max_isolation(:));
[i, j, k, l] = ind2sub(size(max_isolation), I);
fprintf(['Minimum isolation of %.4f:\n' ...
         'Acceleration Gain: %.2f\n' ...
         'Proportional Gain: %.2f\n' ...
         'Dervative Gain: %.2f\n' ...
         'Integral Gain: %.2f\n'], M, x1(i), y1(j), z1(k), a1(l))

save('good_control_info V2', "max_isolation", "acceleration_gains", ...
                                                    "rel_prop_gains", ...
                                                    "rel_deriv_gains", ...
                                                    "rel_int_gains");