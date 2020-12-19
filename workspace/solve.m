%$%
% This method should be inside GraphSLAM class.

function [mu, covar] = solve(obj, tMat, tVec, Mat, Vec)
    % INPUT: tMat, tVec are info mat and vec with tilde
    % with tilde means only the poses, but not the map!
    % OUTPUT: mean value and covariance
    
    % GLOBAL VARIABLES
    numM = size(obj.m, 2);
    
    covar = inv(tMat);
    mu1 = covar * tVec;   % mu 0:t
    
    % compute mean value of map features
    mu2 = zeros(3 * numM);    % (with signature);
    
    for j = 1: numM   
        tau_j = tau_set(j);
        tau_j_idx = find_iMat_idx(obj, tau_j, 'x');   % tau(j): poses
        j_idx = find_iMat_idx(obj, j, 'm');  % j: map
        mu2(3*j-2:3*j) = Mat(j_idx, j_idx) \ (Vec(j_idx) + Mat(j_idx, tau_j_idx) * tVec(tau_j_idx));
    end
    
    mu = [mu1; mu2];
    
end