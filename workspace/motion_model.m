% This function performs the motion model for SLAM target
% Given previous pose mu_{t-1}, control u_t, calculate pose mu_{t}
%
% Inputs:   mu          3X1
%           u           3X1
%
% Output:   mu_bar      3X1

function [mu_bar] = motion_model(mu, u)

    global L % the distance between wheel axles
    alpha = u(1); % steering angle
    v = u(2);    % speed
    delta_t = u(3); % moving time
    fi = mu(3); % heading angle

    mu_bar = mu + 0.001 * delta_t * [v*cos(fi); v*sin(fi);v/L*tan(alpha)];
end