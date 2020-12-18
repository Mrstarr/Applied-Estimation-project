% == j
% i really need interface for getting numX and numM

function [tMat, tVec] = reduce(Mat, Vec)
    % input: info mat and vec
    % =======
    % output: tMat, tVec are info mat and vec with tilde
    % with tilde means only the poses, but not the map!
    
    tMat = Mat;
    tVec = Vec;
    num_f = get_feature_num();  % num of features !!!
    
    for j = 1: num_f
        tau_j = tau_set(j);
        tau_j_idx = find_iMat_idx(tau_j, 'x');   % tau(j): poses
        j_idx = find_iMat_idx(j, 'm');   % m_j idx in info
        all_idx = [tau_j_idx; j_idx];
        
        % update looks strange
        tVec(tau_j_idx) = tVec(tau_j_idx) - tMat(tau_j_idx, j_idx) ...
            / tMat(j_idx, j_idx) * Vec(j_idx);  % (3 * |tau(j)|, 1)
        tMat(tau_j_idx, tau_j_idx) = tMat(tau_j_idx, tau_j_idx) - ...
            tMat(tau_j_idx, j_idx) / tMat(j_idx, j_idx) ...
            * tMat(j_idx, tau_j_idx);
        
    end
    
    global numX
    % remove rows and cols of features
    tMat = tMat(1: 3 * numX, 1: 3 * numX);
    tVec = tVec(1: 3 * numX);
    
end