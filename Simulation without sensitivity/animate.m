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

    theta2_fun = matlabFunction(theta, "Vars", {t});
    theta2_fun = @(x) 0*x;
    theta2 = theta2_fun(t_vec);
    
    x1 = xB + a.l1 * cos(y(1,:) + theta2);
    z1 = zB - a.l1 * sin(y(1,:) + theta2);
    
    %constructing the robotic arm
    h1 = plot(xB(1),zB(1),'ro','MarkerFaceColor','r',"LineWidth",2);
    axis equal
    hold on;
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

    t = get(gca,'Title');
    set(t,'String',strcat('The Robotic Arm Animation (t = ', num2str(floor(t_vec(1))),' seconds)'));
    %title("The Robotic Arm Animation")

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
        set(t,'String',strcat("The Robotic Arm Animation (t = ", num2str(floor(t_vec(ii))), " seconds)"));

        drawnow;
        pause(t_vec(ii)-t_vec(ii-1))
    end
end

