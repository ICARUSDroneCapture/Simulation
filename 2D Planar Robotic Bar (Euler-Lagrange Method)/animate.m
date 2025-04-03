function animate(t_vec,y,d,theta)
    % t_vec: integration time vector [vector]
    % y: integration output vector (solution) [vector]
    % d: deck movement [symbolic function vector]
    % theta: deck pitch angle [symbolic function scalar -- rads]

    %calling the robotic arm constants
    constants;

    pause(1);
    
    %animate
    figure()
    %base
    syms t
    xB = d(1);
    zB = d(3);

    xB = double(subs(xB,t,t_vec));
    zB = double(subs(zB,t,t_vec));
    
    theta2 = double(subs(theta,t,t_vec));

    %link 1
    x1 = xB+platform.l1*cos(y(1,:)+theta2);
    z1 = zB-platform.l1*sin(y(1,:)+theta2);
    
    %constructing the robotic arm
    h1 = plot(xB(1),zB(1),'ro','MarkerFaceColor','r',"LineWidth",2);
    axis equal
    hold on;
    h2 = line([xB(1) x1(1)],[xB(1) z1(1)],'color','b','LineWidth',2);
    h3 = plot(x1(1),z1(1),'go','MarkerFaceColor','g','LineWidth',2);
    axis([-platform.l1-0.1 platform.l1+0.1 -platform.l1-0.1 platform.l1+0.1]);
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

