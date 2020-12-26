% This function plot part of Ground truth of victoria park
% With raw data rotating clockwise theta_0 and translating to (0,0)
% Proportion depends on input percentage [0,1] 
function plot_G(percentage)
    load('aa3_gpsx.mat');
    N = size(Lo_m,1);
    N_p = round(percentage*N);
    theta = -35.5;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    G_x = Lo_m - Lo_m(1);
    G_y = La_m - La_m(1);
    Gt = [G_x G_y]';
    Gt = R * Gt;
    hold on
    plot(Gt(1,1:N_p),Gt(2,1:N_p),'r.','MarkerSize',2);
    hold off
end