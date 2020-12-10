% This function treats motion data to control input
% Input:
%   DR  Motion struct with fields:
%           steering: [61945×1 double] in m/secs.
%           speed: [61945×1 double] in radians.
%           time: [61945×1 double] in millisecs.
% Output: U    3X61944 
function U = motion_dataprocess(DR)
    I = size(DR.time,1);
    U = zeros(3, I-1);
    
    U(1,:) = DR.steering(1:I-1);
    U(2,:) = DR.speed(1:I-1);
    for i = 1:I-1
        U(3,i) = DR.time(i+1)-DR.time(i);
    end
end