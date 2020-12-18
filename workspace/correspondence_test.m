% ==j

function prob = correspondence_test (iMat, iVec, mu, covar, j, k)
    % get mean vector mu and 
    % the path covariance covar(0~t)
    % from GraphSLAM_solve
    % === return ===
    % the posterior prob of mj & mk are the same
    
    jk = find_iMat_idx([j, k], 'm');  % should be implemented in the graph class
    tau_jk_num = tau_set(j, k);  % give a list
%     tau_jk = zeros(3 * size(tau_jk_num(:)), 1);
%     for t = 1: size(tau_jk_num(:))
%         tau_jk = find_iMat_idx(tau_jk_num(t), 'x');
%     end
    tau_jk = find_iMat_idx(tau_jk_num, 'x');  % return the indices
    
    % find tau idx should be implemented in the graph class
    % the set of poses τ (j, k) at which the robot observed feature j and k
    
    iMat_jk = iMat(jk, jk);
    iMat_jkt = iMat(jk, tau_jk);
    iMat_tjk = iMat(tau_jk, jk);
    cov_jk = covar(tau_jk, tau_jk);
    
    % marginalized info
    iMat_margin = iMat_jk - iMat_jkt * cov_jk * iMat_tjk;
    iVec_margin = iMat_margin * mu(jk) ;
    
    % info of the difference variable
    iMat_diff = [1, -1] * iMat_margin * [1; -1];
    iMat_diff_inv = inv(iMat_diff);
    iVec_diff = [1, -1] * iVec_margin;
    mu_diff = iMat_diff_inv * iVec_diff;
    
    prob = 1/sqrt(2*pi*det(iMat_diff_inv)) * ...
        exp(-0.5 * mu_diff' * iMat_diff_inv * mu_diff);
    
end