close all; clc; clear;

simulationParameters;

beta_min = 0.5;
beta_max = 2;
n_max = 15;
isolation = zeros(1,n_max);
max_wave = 2*pi/Tmax;
betas = linspace(beta_min, beta_max, n_max);
beta_wave = 1.46;
[~, wave_idx] = min(abs(betas-beta_wave));

% Computation time: 12.8444 hrs
% acceleration_gains = 10000:1000:25000;
% rel_prop_gains = 0:250:2000;
% rel_deriv_gains = 3000:100:4500;
% rel_int_gains = 0:5:50;
% B_scale = 0:0.2:1;
% r_g = 0.2:0.1:0.5; 
% n_vals = 1:4;
acceleration_gains = 24000;
rel_prop_gains = 250;
rel_deriv_gains = 3000;
rel_int_gains = 10;
B_scale = 1;
r_g = 0.4; 
n_vals = 4;

originalState = warning;

totalIter = length(n_vals) * length(B_scale) * length(r_g) * ...
            length(rel_int_gains) * length(rel_deriv_gains) * ...
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

    warnId = 'MATLAB:ode45:IntegrationTolNotMet';
    warning('error', warnId);

    [a_idx, p_idx, d_idx, i_idx, rg_idx, B_idx, n_idx] = ind2sub([length(acceleration_gains), ...
                                                          length(rel_prop_gains), ...
                                                          length(rel_deriv_gains), ...
                                                          length(rel_int_gains), ...
                                                          length(r_g), ...
                                                          length(B_scale), ...
                                                          length(n_vals)], idx_all);

    isolation = NaN(1, n_max);
    atemp = a;
    atemp.ka = acceleration_gains(a_idx);
    atemp.kp_c = rel_prop_gains(p_idx);
    atemp.kd_c = rel_deriv_gains(d_idx);
    atemp.ki_c = rel_int_gains(i_idx);
    atemp.r_g = r_g(rg_idx);
    atemp.kp_b = B_scale(B_idx) * 3000;
    atemp.kd_b = B_scale(B_idx) * 500;
    atemp.ki_b = B_scale(B_idx) * 200;

    n = n_vals(n_idx);
    atemp.B = @(x) (-((atemp.h_k - 1) / (1 - (atemp.pr_d + atemp.r_k))^n) * ...
                (sign(x - atemp.pr_d) .* (x - (atemp.pr_d + sign(x - atemp.pr_d) * atemp.r_k))).^n + atemp.h_k) .* ...
               (abs(x - atemp.pr_d) > atemp.r_k) + atemp.h_k * (abs(x - atemp.pr_d) <= atemp.r_k);

    atemp.C = @(x) ((-1 / (1 - (atemp.pr_d + atemp.r_g))^n) * ...
                (sign(x - atemp.pr_d) .* (x - (atemp.pr_d + sign(x - atemp.pr_d) * atemp.r_g))).^n + 1) .* ...
               (abs(x - atemp.pr_d) > atemp.r_g) + 1 * (abs(x - atemp.pr_d) <= atemp.r_g);

    % Deck wave angular frequencies to test
    betas = linspace(beta_min, beta_max, n_max);

    for idx = 1:n_max
        beta = betas(idx);
        atemp.d = @(t) alpha * cos(beta * t) + hdeck;
        atemp.d_dot = @(t) -beta * alpha * sin(beta * t);
        atemp.d_ddot = @(t) -beta^2 * alpha * cos(beta * t);

        [t, s] = controlSystemsTesting(atemp);

        if plotSim & (idx == wave_idx)
            t_wave_temp{idx_all} = t
            s_wave_temp{idx_all} = s
            a_wave_temp{idx_all} = atemp;
        end


        t_test0 = 70;
        if t(end) == atemp.tspan(end)
            x_I = s(t > t_test0, 1);
            x_D = atemp.d(t(t > t_test0));
            isolation(idx) = calculateAverageIsolation(x_I, x_D);
            if isolation(idx) < 0
                isolation(idx) = NaN;
            end
        end
    end

    if plotSim
        isolation_wave_temp{idx_all} = isolation;
    end

    if sum(isnan(isolation)) == 0
        max_isolation_temp{idx_all} = max(isolation(betas >= max_wave));
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
    plotResponse(t_wave, s_wave, betas, isolation_wave, a_wave, max_isolation)
end

% Reconstruct max_isolation after parfor
max_isolation = NaN(length(acceleration_gains), length(rel_prop_gains), ...
                    length(rel_deriv_gains), length(r_g), length(B_scale), length(n_vals));
for idx_all = 1:totalIter
    [a_idx, p_idx, d_idx, i_idx, rg_idx, B_idx, n_idx] = ind2sub([length(acceleration_gains), ...
                                                          length(rel_prop_gains), ...
                                                          length(rel_deriv_gains), ...
                                                          length(rel_int_gains), ...
                                                          length(r_g), ...
                                                          length(B_scale), ...
                                                          length(n_vals)], idx_all);
    max_isolation(a_idx, p_idx, d_idx, i_idx, rg_idx, B_idx, n_idx) = max_isolation_temp{idx_all};
end

toc

warning(originalState);

save('good_control_info_integral_V5', "max_isolation", "acceleration_gains", ...
                                                    "rel_prop_gains", ...
                                                    "rel_deriv_gains", ...
                                                    "rel_int_gains", ...
                                                    "r_g", ...
                                                    "B_scale", ...
                                                    "n_vals");
