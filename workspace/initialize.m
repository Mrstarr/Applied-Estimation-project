% This function performs GRAPH SLAM initialization 
% It takes the controls u1:t as input and outputs sequence of pose estimates ?0:t
% ?0:t constructs the graph node
%
% Input:    u           3XN
% Output:   mu_full     3X(N+1)
function mu_full = initialize(u)
    I = size(u,2) % size of motion nodes
    mu_full = zeros(3,I+1);
    mu_full(:,1) = zeros(3,1);
    
    for i = 1:I
        mu_full(:,i+1) = motion_model(mu_full(:,i),u(:,i));
    end
end