function showgraph()
    x = load("state.mat","-mat");
    m = load("landmark.mat","-mat");
    hold on
%     x0 = [-67.6493; -41.7142; 35.5*pi/180];
    plot(x.state(1,:),x.state(2,:),'k');
    plot(m.landmark(1,:),m.landmark(2,:),'b.','MarkerSize', 5);
end