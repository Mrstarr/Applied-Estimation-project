function [ G ]= calculate_jacobian(x,u)
    global vehicle
    a = vehicle.a;
    b = vehicle.b;
    L = vehicle.L;
    
    v = u(1);    % speed
    alpha = u(2); % steering angle
    T = u(3); % moving time
    fi = x(3); % heading angle
    
    G = [1 0 -T*v*(sin(fi) + 1/L*tan(alpha)*(a*cos(fi) - b*sin(fi))); 
         0 1  T*v*(cos(fi) - 1/L*tan(alpha)*(a*sin(fi) + b*cos(fi)));
         0 0  1];
end
