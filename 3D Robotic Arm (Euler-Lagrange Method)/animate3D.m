function animate3D(t_vec,y,platform,d,theta_D)
    
    %calling the robotic arm constants
    % constants;

    pause(1);
    
    %animate
    figure()

    %base translations
    syms t
    xB = double(subs(d(1),t,t_vec));
    yB = double(subs(d(2),t,t_vec));
    zB = double(subs(d(3),t,t_vec));

    %base rotations
    d1 = double(subs(theta_D(1),t,t_vec));
    d2 = double(subs(theta_D(2),t,t_vec));
    d3 = double(subs(theta_D(3),t,t_vec));

    %link 1
    x1 = xB + platform.l1*cos(d2).*cos(d3).*cos(y(1,:)).*cos(y(3,:)) - platform.l1*cos(y(1,:)).*sin(y(3,:)).*(cos(d1).*sin(d3) - cos(d3).*sin(d1).*sin(d2)) - platform.l1*sin(y(1,:)).*(sin(d1).*sin(d3) + cos(d1).*cos(d3).*sin(d2));
    y1 = yB + platform.l1*sin(y(1,:)).*(cos(d3).*sin(d1) - cos(d1).*sin(d2).*sin(d3)) + platform.l1*cos(y(1,:)).*sin(y(3,:)).*(cos(d1).*cos(d3) + sin(d1).*sin(d2).*sin(d3)) + platform.l1*cos(d2).*cos(y(1,:)).*cos(y(3,:)).*sin(d3);
    z1 = zB + platform.l1*cos(d2).*cos(y(1,:)).*sin(d1).*sin(y(3,:)) - platform.l1*cos(y(1,:)).*cos(y(3,:)).*sin(d2) - platform.l1*cos(d1).*cos(d2).*sin(y(1,:));

    %link 2
    x2 = xB + cos(d2).*cos(d3).*(platform.l1*cos(y(1,:)).*cos(y(3,:)) + platform.l2*cos(y(3,:)).*cos(y(1,:) + y(2,:))) - (sin(d1).*sin(d3) + cos(d1).*cos(d3).*sin(d2)).*(platform.l1*sin(y(1,:)) + platform.l2*sin(y(1,:) + y(2,:))) - (cos(d1).*sin(d3) - cos(d3).*sin(d1).*sin(d2)).*(platform.l1*cos(y(1,:)).*sin(y(3,:)) + platform.l2*sin(y(3,:)).*cos(y(1,:) + y(2,:)));
    y2 = yB + (cos(d1).*cos(d3) + sin(d1).*sin(d2).*sin(d3)).*(platform.l1*cos(y(1,:)).*sin(y(3,:)) + platform.l2*sin(y(3,:)).*cos(y(1,:) + y(2,:))) + (cos(d3).*sin(d1) - cos(d1).*sin(d2).*sin(d3)).*(platform.l1*sin(y(1,:)) + platform.l2*sin(y(1,:) + y(2,:))) + cos(d2).*sin(d3).*(platform.l1*cos(y(1,:)).*cos(y(3,:)) + platform.l2*cos(y(3,:)).*cos(y(1,:) + y(2,:)));
    z2 = zB + cos(d2).*sin(d1).*(platform.l1*cos(y(1,:)).*sin(y(3,:)) + platform.l2*sin(y(3,:)).*cos(y(1,:) + y(2,:))) - sin(d2).*(platform.l1*cos(y(1,:)).*cos(y(3,:)) + platform.l2*cos(y(3,:)).*cos(y(1,:) + y(2,:))) - cos(d1).*cos(d2).*(platform.l1*sin(y(1,:)) + platform.l2*sin(y(1,:) + y(2,:))); 

    %constructing the robotic arm
    axis equal
    h1 = plot3(xB(1),yB(1),zB(1),'ro','MarkerFaceColor','r',"LineWidth",2);
    hold on;
    h2 = line([xB(1) x1(1)],[yB(1) y1(1)],[zB(1) z1(1)],'color','b','LineWidth',2);
    h3 = plot3(x1(1),y1(1),z1(1),'go','MarkerFaceColor','g','LineWidth',2);
    h4 = line([x1(1) x2(1)],[y1(1) y2(1)],[z1(1) z2(1)],'color','k','LineWidth',2);
    h5 = plot3(x2(1),y2(1),z2(1),'m^','MarkerFaceColor','m','LineWidth',2);
    axis([-platform.l1-platform.l2-0.02-6 platform.l1+platform.l2+0.02+6 -platform.l1-platform.l2-0.02 platform.l1+platform.l2+0.02 -platform.l1-platform.l2-0.02-0.5 platform.l1+platform.l2+0.02+0.5]);
    xlabel('x-axis [m]');
    ylabel('y-axis [m]');
    zlabel('z-axis [m]');

    %plotting the xy surface
    x = xB(1)-(platform.l1+platform.l2+0.02):0.3:xB(1)+(platform.l1+platform.l2+0.02);
    y = (yB(1)-(platform.l1+platform.l2+0.02)):0.3:(yB(1)+(platform.l1+platform.l2+0.02));
    a = sin(d1(1))*sin(d3(1))+cos(d1(1))*cos(d3(1))*sin(d2(1));
    b = cos(d1(1))*sin(d2(1))*sin(d3(1))-cos(d3(1))*sin(d1(1));
    c = cos(d1(1))*cos(d2(1));
    [X,Y] = meshgrid(x,y);
    Z = ((a*(xB(1)-X)+b*(yB(1)-Y))/(c)+zB(1));

    h6 = surf(X,Y,Z);
    h6.FaceColor = [0.6350 0.0780 0.1840];

    %plotting the inertial surface at z = 0
    xI = -10:0.3:10;
    yI = -10:0.3:10;
    [XI, YI] = meshgrid(xI,yI);
    ZI = zeros(size(XI));
    
    h7 = surf(XI,YI,ZI);
    h7.FaceAlpha = 0.1;
    h7.EdgeAlpha = 0.1;
    
    t = get(gca,'Title');
    set(t,'String',strcat('The Robotic Arm Animation (t = ', num2str(floor(t_vec(1))),' seconds)'));

    %title("The Robotic Arm Animation")
    view(-10,10)
    grid on
    for ii=2:1:length(x1)
        set(h1,'XData',xB(ii));
        set(h1,'YData',yB(ii));
        set(h1,'ZData',zB(ii));
        set(h2,'XData',[xB(ii) x1(ii)]);
        set(h2,'YData',[yB(ii) y1(ii)]);
        set(h2,'ZData',[zB(ii) z1(ii)]);
        set(h3,'XData',x1(ii));
        set(h3,'YData',y1(ii));
        set(h3,'ZData',z1(ii));
        set(h4,'XData',[x1(ii) x2(ii)]);
        set(h4,'YData',[y1(ii) y2(ii)]);
        set(h4,'ZData',[z1(ii) z2(ii)]);
        set(h5,'XData',x2(ii));
        set(h5,'YData',y2(ii));
        set(h5,'ZData',z2(ii));
        x = xB(ii)-(platform.l1+platform.l2+0.02):0.3:xB(ii)+(platform.l1+platform.l2+0.02);
        y = (yB(ii)-(platform.l1+platform.l2+0.02)):0.3:(yB(ii)+(platform.l1+platform.l2+0.02));
        [X,Y] = meshgrid(x,y);
        a = sin(d1(ii))*sin(d3(ii))+cos(d1(ii))*cos(d3(ii))*sin(d2(ii));
        b = cos(d1(ii))*sin(d2(ii))*sin(d3(ii))-cos(d3(ii))*sin(d1(ii));
        c = cos(d1(ii))*cos(d2(ii));
        Z = ((a*(xB(ii)-X)+b*(yB(ii)-Y))/(c)+zB(ii));

        set(h6,'XData',X);
        set(h6,'YData',Y);
        set(h6,'ZData',Z);

        set(t,'String',strcat("The Robotic Arm Animation (t = ", num2str(floor(t_vec(ii))), " seconds)"));

        drawnow;
        % pause(t_vec(ii)-t_vec(ii-1))
    end
    hold off;
end

