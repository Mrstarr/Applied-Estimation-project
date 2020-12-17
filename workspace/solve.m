% == j
% what's the size of miu?

function [mu, covar] = solve(tMat, tVec, Mat, Vec)
    % tMat, tVec are info mat and vec with tilde
    % with tilde means only the poses, but not the map!
    % === return ===
    % mean value and covariance
    
    covar = inv(tMat);
    mu1 = covar * tVec;   % mu 0:t
    num_f = get_feature_num();   % =========!!!!!!!! should get length
    
    % compute mean value of map features
    mu2 = zeros(3 * num_f);    % (with signature);
    
    for j = 1: length(map)   
        tau_j = tau_set(j);
        tau_j_idx = find_iMat_idx(tau_j, 'x');   % tau(j): poses
        j_idx = find_iMat_idx(j, 'm');  % j: map
        mu2(3*j-2:3*j) = Mat(j_idx, j_idx) \ (Vec(j_idx) + Mat(j_idx, tau_j_idx) * tVec(tau_j_idx));
    end
    
    mu = [mu1; mu2];
    
end