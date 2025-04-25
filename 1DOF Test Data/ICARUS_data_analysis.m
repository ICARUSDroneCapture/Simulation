clear
close all
clc

files = dir(fullfile("ICARUSDATA","*.csv"));
titles = ["Fall from 90 Degrees",
            "Acceleration Control, ka = 50",
            "Relative Position Control, kp = 100, starting angle = 0",
            "Relative Position Control, kp = 100, starting angle = 90",
            "Relative Position Control, kp = 100, kd = 7, ki = 1",
            "Relative Position Control, kp = 100, kd = 7, ki = 1, starting angle = 90",
            "Acceleration Control, ka = 50, test #2", 
            "Acceleration Control, ka = 500", ];
simData = [0, 'video2SimData',0,0,'video5SimData',0,0,0];

start_time = [8.35, 11.32, 0, 0, 10.1, 0, 0, 0];

outputFolder = "plots"; 
if ~exist(outputFolder,"dir")
    mkdir(outputFolder);
end

% for k = 1:numel(files)
for k = 1
    fname = fullfile(files(k).folder, files(k).name);
    M     = readmatrix(fname);      % N×3 numeric: [t x y]
    t     = M(:,1);
    xd     = M(:,2) * 0.0254;
    y     = M(:,3) * 0.0254;
    [~, name, ~] = fileparts(files(k).name);

    %% 1) Y vs t
    load(sprintf('video%dSimData',k))
    h1 = figure;
    idx = t > start_time(k);
    t_act = t(idx) - start_time(k);
    y_act = y(idx);
    y_act = y_act - max(y_act);
    plot(t_act, y_act, 'LineWidth', 2)
    hold on
    c = 0.73;
    t_max = t_act(end);
    simPos = zEE*c - max(zEE*c) - 0.005;
    plot(x(x < t_max), simPos(x < t_max), 'LineWidth', 1)
    sim_interp = interp1(x(x < t_max), simPos(x < t_max), t_act, 'linear', 'extrap');
    data_range = max(sim_interp) - min(sim_interp);
    fprintf('Average Residual for test %d: %.3f\n', k, mean(abs(sim_interp - y_act)./data_range))
    title(titles(k), 'Interpreter','none')
    xlabel("Time (s)")
    ylabel("Y Position (m)")
    legend('Test Data', 'Simulation', 'Location','northeast')
    grid on
    saveas(h1, fullfile(outputFolder, name + "_y_vs_t.png"))

    %% 2) X vs t
    h2 = figure;
    plot(t, xd, 'LineWidth', 2)
    title(titles(k), 'Interpreter','none')
    xlabel("Time (s)")
    ylabel("X Position (inches)")
    grid on
    saveas(h2, fullfile(outputFolder, name + "_x_vs_t.png"))

    %% 3) Comet (trajectory)
    h3 = figure;
    comet(xd, y)
    xlabel("X Position (inches)")
    ylabel("Y Position (inches)")
    title(titles(k), 'Interpreter','none')
    grid on
    drawnow   % ensure the final frame is rendered
    % snapshot the axes to a PNG
    saveas(h3, fullfile(outputFolder, name + "_trajectory.png"))

    % (optional) close them to keep your desktop clean
    % close([h1,h2,h3])
end
