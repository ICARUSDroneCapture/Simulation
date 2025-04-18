function plot2_test(t,y)
% plot2 plots one figure of q3, Dq3, q1, Dq1, q2, and Dq2
%
% Inputs:
%   t           - time
%   y           - outputs
%
% Outputs: one figure


    figure()
    subplot(3,2,1)
    plot(t,y(3,:))
    xlabel("Time [s]")
    ylabel("q1 [rad]")
    title("q1 Over Time")
    yline(0,'r')
    yline(-pi,'r')
    
    subplot(3,2,2)
    plot(t,y(4,:))
    xlabel("Time [s]")
    ylabel("Dq1 [rad/s]")
    title("Dq1 Over Time")   

    subplot(3,2,3)
    plot(t,y(1,:))
    xlabel("Time [s]")
    ylabel("q2 [rad]")
    title("q2 Over Time")
    yline(pi,'r')
    yline(-pi,'r')
    
    subplot(3,2,4)
    plot(t,y(2,:))
    xlabel("Time [s]")
    ylabel("Dq2 [rad/s]")
    title("Dq2 Over Time") 

    subplot(3,2,5)
    plot(t,y(5,:))
    xlabel("Time [s]")
    ylabel("q3 [rad]")
    title("q3 Over Time")
    yline(pi,'r')
    yline(-pi,'r')
    
    subplot(3,2,6)
    plot(t,y(6,:))
    xlabel("Time [s]")
    ylabel("Dq3 [rad/s]")
    title("Dq3 Over Time")
end

