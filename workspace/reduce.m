%$%
% This method should be inside GraphSLAM class.

function [tMat, tVec] = reduce(obj, Mat, Vec)
    % input: info mat and vec
    % =======
    % output: tMat, tVec are info mat and vec with tilde
    % with tilde means only the poses, but not the map!
    
    % GLOBAL VARIABLES
    numX = size(obj.x, 2);
    numM = size(obj.m, 2);

    % init
    tMat = Mat;
    tVec = Vec;
    
    for j = 1: numM
        tau_j = tau_set(j);
        tau_j_idx = find_iMat_idx(obj, tau_j, 'x');   % tau(j): poses
        j_idx = find_iMat_idx(obj, j, 'm');   % m_j idx in info
        
        % update looks strange
        tVec(tau_j_idx) = tVec(tau_j_idx) - tMat(tau_j_idx, j_idx) ...
            / tMat(j_idx, j_idx) * Vec(j_idx);  % (3 * |tau(j)|, 1)
        tMat(tau_j_idx, tau_j_idx) = tMat(tau_j_idx, tau_j_idx) - ...
            tMat(tau_j_idx, j_idx) / tMat(j_idx, j_idx) ...
            * tMat(j_idx, tau_j_idx);
        
    end
    
    % remove rows and cols of features
    tMat = tMat(1: 3 * numX, 1: 3 * numX);
    tVec = tVec(1: 3 * numX);
    
end