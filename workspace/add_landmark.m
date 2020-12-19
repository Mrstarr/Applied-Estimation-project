% This function takes the moment when measurement zt is received
% Then initialize the features in the map
function [ldm] = add_landmark(x,z)
    if isempty(z)
        ldm = [];
    else
        ldm = x(1:2) + [z(1,:).*cos(z(2,:) + x(3) - pi/2); z(1,:).*sin(z(2,:) + x(3) - pi/2)];
        ldm = [ldm;z(3,:)];
    end
end