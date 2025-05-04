function a_I = Rotate_S_I(a_S, theta, phi, psi)
% Rotate_I_S  Rotates acceleration vector from sensor frame to inertial frame.
%   a_I = Rotate_I_S(a_S, theta, phi, psi) converts a_S vector to a_I vector based on all angles.
%
% Inputs:
%   a_S     : a (column) vector representing the sensor frame acceleration
%             vector in sensor x, y, and x. Eg. [a_Sx; a_Sy; a_Sz] at some point in time.
%
%   theta   : (elevation) angle of platform, aka platform pitch, at some
%             point in time
%
%   phi     : (bank) angle of platform, aka platform roll, at some point in
%             time
%
%   psi     : (azmuth) angle of platform, aka platform yaw, at some point
%             in time
%
% Outputs:
%   a_I     : a (column) vector representing the inertial acceleration vector 
%             in a, y, and z. Eg. [a_Ix; a_Iy; a_Iz] at some point in time.

    s = @(x) sin(x);
    c = @(x) cos(x);

    R_I_1 = [c(psi)  s(psi) 0;
         -s(psi) c(psi) 0;
         0       0      1];

    R_1_2 = [1      0       0;
             0      c(phi)  s(phi);
             0      -s(phi) c(phi)];

    R_2_3 = [c(theta)  0    s(theta);
             0         1    0;
             -s(theta) 0    c(theta)];

    % R_I_1 = [c(psi)     s(psi)      0;
    %          -s(psi)    c(psi)      0;
    %          0          0           1];
    % 
    % R_1_2 = [c(theta)   0           -s(theta);
    %          0          1           s(phi);
    %          s(theta)   0           c(theta)];
    % 
    % R_2_3 = [1          0           0;
    %          0          c(phi)      s(phi);
    %          0          -s(phi)     c(phi)];

    
    R_I_S = R_2_3 * R_1_2 * R_I_1;

    a_I = R_I_S' * a_S;

end