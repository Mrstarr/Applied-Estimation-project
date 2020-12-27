%%
% NOTE: The path of whole project directory should be added.
% And its sub directories.

%% Global variables
global vehicle x0 noise I 
vehicle.L = 2.83;
vehicle.a = 0.95;
vehicle.b = 0.5;
vehicle.H = 0.76;
x0 = [-67.6493; -41.7142; 35.5*pi/180];

noise.R = diag([0.05 0.05 0.001]); % (x, y, th)  [0.05 0.05 0.001]
noise.Q = diag([1 0.01 0.1]);   % (range, angle, signature)     [1 0.01]
I = 2500;

%% test: DELETE outliers
global DEL_THRESH
DEL_THRESH = 2;   % thresh <= 1 means we don't delete
% or it means that any idx below this thresh will be deleted

%%
% read dead reckoning data and laser data
load('VictoriaParkSinSincronizar.mat')
gs = GraphSLAM();
gs.initialize(time,TLsr,u,zt);
gs.plotgraph("r");
gs.run(3);
