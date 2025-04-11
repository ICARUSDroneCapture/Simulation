function plot1(t,y,ref_q1)
% plot1 plots one figure of q1, Dq1, q2, and Dq2 for positional control
%
% Inputs:
%   t           - time
%   y           - outputs
%   ref_q1      - reference q1 for positional control
%   ref_q2      - reference q2 for positional control
%
% Outputs: one figure


    figure()
    subplot(1,2,1)
    plot(t,y(1,:))
    xlabel("Time [s]")
    ylabel("q1 [rad]")
    title("q1 Over Time")
    yline(pi,'r')
    yline(-pi,'r')
    yline(ref_q1, 'b--')
    
    subplot(1,2,2)
    plot(t,y(1,:))
    xlabel("Time [s]")
    ylabel("Dq1 [rad/s]")
    title("Dq1 Over Time")
    
end

