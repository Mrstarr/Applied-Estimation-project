% Global variables
global L 
L = 2.83


%%
DR = load('aa3_dr.mat','-mat');
u = motion_dataprocess(DR);
 mu_full = initialize(u);
plot(mu_full(1,:),mu_full(2,:))
