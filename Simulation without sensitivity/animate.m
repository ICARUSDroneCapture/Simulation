function animate(t_vec,y,d,theta, a)
    % t_vec: integration time vector [vector]
    % y: integration output vector (solution) [vector]
    % d: deck movement [symbolic function vector]
    % theta: deck pitch angle [symbolic function scalar -- rads]
    % a: structure of simulation parameters

    pause(1);
    
    %animate
    figure()
    %base
    syms t
    
    d_fun = matlabFunction(d, "Vars", {t});
    
    B = zeros(3, length(t_vec));
    for i = 1:length(t_vec)
        B(:, i) = d_fun(t_vec(i));
    end
    xB = B(1, :);
    zB = B(3, :);

    theta2 = theta(t_vec);
    
    x1 = xB + a.l1 * cos(y(1,:) + theta2);
    z1 = zB - a.l1 * sin(y(1,:) + theta2);

    x_pr = @(x, q) a.l1 * cos(x+q);
    z_pr = @(x, q) -a.l1 * sin(x+q);

    b_range_n = linspace(-17*pi/32, -pi/2, 5);
    bl_range_n = linspace(-pi/2, -pi/4, 25);
    c_range = linspace(-pi/4, pi/4, 50);
    bl_range_p = linspace(pi/4, pi/2, 25);
    b_range_p = linspace(pi/2, 17*pi/32, 5);

    theta_a = theta2(1) + a.q1_ref;
    h8 = fill([0 x_pr(theta_a, b_range_n) 0], ...
              [0 z_pr(theta_a, b_range_n) 0], "r", ...
              'FaceColor', "#ea9998", 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    h9 = fill([0 x_pr(theta_a, bl_range_n) 0], ...
              [0 z_pr(theta_a, bl_range_n) 0], "r", ...
              'FaceColor', "#ffe59a", 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    h10 = fill([0 x_pr(theta_a, c_range) 0], ...
              [0 z_pr(theta_a, c_range) 0], "r", ...
              'FaceColor', "#b6d7a8", 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    h11 = fill([0 x_pr(theta_a, bl_range_p) 0], ...
              [0 z_pr(theta_a, bl_range_p) 0], "r", ...
              'FaceColor', "#ffe59a", 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    h12 = fill([0 x_pr(theta_a, b_range_p) 0], ...
              [0 z_pr(theta_a, b_range_p) 0], "r", ...
              'FaceColor', "#ea9998", 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    %constructing the robotic arm
    h1 = plot(xB(1),zB(1),'ro','MarkerFaceColor','r',"LineWidth",2);
    axis equal
    h2 = line([xB(1) x1(1)],[xB(1) z1(1)],'color','b','LineWidth',2);
    h3 = plot(x1(1),z1(1),'go','MarkerFaceColor','g','LineWidth',2);
    axis([-a.l1-0.1 a.l1+0.1 -a.l1-0.1 a.l1+0.1]);
    xlabel('x-axis [m]');
    ylabel('z-axis [m]');

    yline(0,'k--','Label','Inertial')

    xspan = [xB(1)-0.5*(sqrt(4/(1+(tan(theta2(1)))^2))) 0];
    deckLine = [tan(theta2(1))*(xspan(1)-xB(1))+zB(1) tan(theta2(1))*(xspan(2)-xB(1))+zB(1)];
    h4 = line(xspan,deckLine,'LineWidth',2);
    h4.Color = [0.6350 0.0780 0.1840];

    q_d = a.q1_ref;
    q_max = a.q1_ref + a.w/2;
    q_min = a.q1_ref - a.w/2;
    h5 = line([xB(1) x_pr(theta2(1),q_d)], ...
              [zB(1) z_pr(theta2(1),q_d)], 'LineStyle', '--', ...
                                            'LineWidth', 2);
    h6 = line([xB(1) x_pr(theta2(1),q_max)], ...
              [zB(1) z_pr(theta2(1),q_max)], 'LineStyle', '--', ...
                                             'Color', 'k', ...
                                             'LineWidth', 2);
    h7 = line([xB(1) x_pr(theta2(1),q_min)], ...
              [zB(1) z_pr(theta2(1),q_min)], 'LineStyle', '--', ...
                                             'Color', 'k', ...
                                             'LineWidth', 2);

    t = get(gca,'Title');
    set(t,'String',strcat('The Robotic Arm Animation (t = ', num2str(floor(t_vec(1))),' seconds)'));
    %title("The Robotic Arm Animation")

    legend([h2, h4, h5, h6], ...
        {'1-DOF Arm','Deck','Desired Position', 'Control Boundaries'}, ...
        'FontSize', 8, 'Location','southwest')

    for ii=2:1:length(x1)
        set(h1,'XData',xB(ii));
        set(h1,'YData',zB(ii));
        set(h2,'XData',[xB(ii) x1(ii)]);
        set(h2,'YData',[zB(ii) z1(ii)]);
        set(h3,'XData',x1(ii));
        set(h3,'YData',z1(ii));
        xspan = [xB(ii)-0.5*(sqrt(4/(1+(tan(theta2(ii)))^2))) 0];
        set(h4,'XData',xspan);
        set(h4,'YData',[-tan(theta2(ii))*(xspan(1)-xB(ii))+zB(ii) -tan(theta2(ii))*(xspan(2)-xB(ii))+zB(ii)]);
        set(h5,'XData',[xB(ii) x_pr(theta2(ii),q_d)]);
        set(h5,'YData',[zB(ii) z_pr(theta2(ii),q_d)]);
        set(h6,'XData',[xB(ii) x_pr(theta2(ii),q_max)]);
        set(h6,'YData',[zB(ii) z_pr(theta2(ii),q_max)]);
        set(h7,'XData',[xB(ii) x_pr(theta2(ii),q_min)]);
        set(h7,'YData',[zB(ii) z_pr(theta2(ii),q_min)]);
        theta_a = theta2(ii) + a.q1_ref;
        set(h8,'XData',[0 x_pr(theta_a, b_range_n) 0]);
        set(h8,'YData',[0 z_pr(theta_a, b_range_n) 0]);
        set(h9,'XData',[0 x_pr(theta_a, bl_range_n) 0]);
        set(h9,'YData',[0 z_pr(theta_a, bl_range_n) 0]);
        set(h10,'XData',[0 x_pr(theta_a, c_range) 0]);
        set(h10,'YData',[0 z_pr(theta_a, c_range) 0]);
        set(h11,'XData',[0 x_pr(theta_a, bl_range_p) 0]);
        set(h11,'YData',[0 z_pr(theta_a, bl_range_p) 0]);
        set(h12,'XData',[0 x_pr(theta_a, b_range_p) 0]);
        set(h12,'YData',[0 z_pr(theta_a, b_range_p) 0]);
        set(t,'String',strcat("The Robotic Arm Animation (t = ", num2str(floor(t_vec(ii))), " seconds)"));

        drawnow;
        % pause(t_vec(ii)-t_vec(ii-1))
    end
end

