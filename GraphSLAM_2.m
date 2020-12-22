classdef GraphSLAM < handle
    %POSEGRAPH A class for doing pose graph optimization
    
    properties (SetAccess = private)
        u    % control 1:t
        x     % motion state 1:t
        m     % maps 
        z     % effective measurement
        H     % Information matrix
        b     % Information vector
        tau % the set of all poses xt at which landmark was observed
        n_node  % Number of nodes in graph
        n_edge  % Number of edges in graph
        n_ldm   % Number of landmarks
    end  % properties set private
    
    properties (Dependent = true)
        
        % pose    % Poses of all nodes
    end  % properties dependent
    
    methods

        % initialize the pose graph
        function initialize(obj,time,TLsr,u,zt)
            I = round(size(u,2)/20);
            ILsr = size(zt,2);
            
            % u & x
            obj.u = zeros(3, I-1);
            obj.x = zeros(3,I);
            % get control data
            obj.u(1,:) = u(1,1:I-1); % motion speed
            obj.u(2,:) = u(2,1:I-1); % steering angle
            obj.u(3,:) = diff(time(1,1:I)); % time difference           
            
            global x0
            obj.x(:,1) = x0;
            
            % tranverse all control and measurement 
            j = 1;
            obj.n_ldm = 0;
            obj.z.m = []; % measurement 
            obj.z.idx = []; % time index i
            obj.z.c = []; % correspondence
            %obj.m =[]; % landmark feature
            %obj.m.tau = []; % set of poses that observes the landmark
            for i = 1:size(obj.u,2)
                obj.x(:,i+1) = predict(obj.x(:,i),obj.u(:,i));
                if j <= ILsr
                    if TLsr(j) < time(i)
                        % add new measurement
                        obj.z.m = [obj.z.m zt{j}];  
                        obj.z.idx = [obj.z.idx i * ones(1,size(zt{j},2))]; 
                        obj.z.c = [obj.z.c (obj.n_ldm + 1):obj.n_ldm + size(zt{j},2)]; 
                       
                        % add new landmark
                        ldm = add_landmark(obj.x(:,i),zt{j});   
                        obj.m = [obj.m ldm];
                        for q = (obj.n_ldm + 1) : (obj.n_ldm + size(zt{j},2))
                            obj.tau{q} = [ i ];
                        end
                        obj.n_ldm = obj.n_ldm + size(zt{j},2);
                        j = j + 1;
                    end
                end     
            end
            % update node and edge number
            fprintf("Updated Node Number = %d\n",obj.n_node)
            fprintf("Pose Graph Initialization Complete\n")
            fprintf("Number of States : %d\n", size(obj.x,2));
            fprintf("Number of Landmarks : %d\n", size(obj.m,2));
            
        end
        
        
        function linearize(obj)
            global noise
            R = noise.R;
            Q = noise.Q;
            % initialize H and b
            obj.n_node = size(obj.x,2) + size(obj.m,2);
            obj.H = zeros(3*obj.n_node,3*obj.n_node);
            obj.b = zeros(obj.n_node*3,1);

            H0 = 1e8 * eye(3);
            obj.H(1:3,1:3) = H0;
            
            % tranverse all control 
            for i = 1:size(obj.u,2)
                % calculate jacobian and H_ii,H_ij,H_jj
                % add u
                G = calculate_jacobian(obj.x(:,i),obj.u(:,i));                     
                H_ii = inv(R);
                H_ij = - R \ G;
                H_jj = G' / R * G ;
                b_i = R \ (obj.x(:,i+1) - G * obj.x(:,i));
                b_j = - G'/ R * (obj.x(:,i+1) - G * obj.x(:,i));

                i_idx = (3*(i-1)+1):3*i;
                j_idx = (3*i+1):3*(i+1);
                
                obj.H(i_idx,i_idx) = obj.H(i_idx,i_idx) + H_ii;
                obj.H(i_idx,j_idx) = obj.H(i_idx,j_idx) + H_ij;
                obj.H(j_idx,i_idx) = obj.H(j_idx,i_idx) + H_ij';
                obj.H(j_idx,j_idx) = obj.H(j_idx,j_idx) + H_jj;
                obj.b(i_idx) = obj.b(i_idx) + b_i;
                obj.b(j_idx) = obj.b(j_idx) + b_j;
            end
           
            % tranverse all measurement 
            for j = 1:size(obj.z.m,2)
                
                cj = obj.z.c(j);
                t = obj.z.idx(j);
                dlt = [obj.m(1,cj) - obj.x(1,t);
                         obj.m(2,cj) - obj.x(2,t)];
                q = dlt'*dlt;
                z_hat = [sqrt(q);atan2(dlt(2),dlt(1))-obj.x(3,t)+pi/2;obj.m(3,cj)];
                H_1 = 1/q * [ -sqrt(q) * dlt(1) -sqrt(q) * dlt(2) 0;
                              dlt(2)  -dlt(1)  -q;
                              0 0 0];
                H_2 = 1/q * [ sqrt(q) * dlt(1) sqrt(q) * dlt(2) 0;
                              -dlt(2)  dlt(1)  0;
                              0 0 q];

                % update information matrix       
                t_idx = (3*(t-1)+1):3*t;
                j_idx = (3*(cj-1)+1+3*size(obj.x,2)):(3*cj+3*size(obj.x,2));
                H_11 = H_1' / Q * H_1;
                H_12 = H_1' / Q * H_2;
                H_21 = H_2' / Q * H_1;
                H_22 = H_2' / Q * H_2;
                b1 = H_1' / Q * (obj.z.m(:,j) - z_hat + [H_1 H_2]*[obj.x(:,t); obj.m(:,cj)]);
                b2 = H_2' / Q * (obj.z.m(:,j) - z_hat + [H_1 H_2]*[obj.x(:,t); obj.m(:,cj)]);

                obj.H(t_idx,t_idx) = obj.H(t_idx,t_idx) + H_11;
                obj.H(t_idx,j_idx) = obj.H(t_idx,j_idx) + H_12;
                obj.H(j_idx,t_idx) = obj.H(j_idx,t_idx) + H_21;
                obj.H(j_idx,j_idx) = obj.H(j_idx,j_idx) + H_22;
                obj.b(t_idx) = obj.b(t_idx) + b1;
                obj.b(j_idx) = obj.b(j_idx) + b2;
            end
            
%             H_sparse = sparse(obj.H);
%             mu = H_sparse \ obj.b;
%             Mu = reshape(mu, 3, size(mu,1)/3);
%             numX = size(obj.x, 2)
%             plot(Mu(1,1:numX),Mu(2,1:numX),'b')
            fprintf("Linearization Complete\n")
            fprintf("Size of Information matrix: %d * %d\n", size(obj.H))
        end
        
        
        function [tMat, tVec] = reduce(obj)
            tMat = obj.H;
            tVec = obj.b;
            numM = size(obj.m, 2);
            
            for j =1:numM
                for g = obj.tau{j}
                    % find corresponding index of j and x
                    idx_x = (3*(g-1)+1):3*g;
                    idx_j = (3*(j-1)+1+3*size(obj.x,2)):(3*j+3*size(obj.x,2));
                    % reduce H and b
                    tVec(idx_x) = tVec(idx_x)- tMat(idx_x,idx_j)/tMat(idx_j,idx_j)*tVec(idx_j);
                    tMat(idx_x,idx_x) = tMat(idx_x,idx_x)- tMat(idx_x,idx_j)/tMat(idx_j,idx_j)*tMat(idx_j,idx_x);
                end
                % remove rows/columns corresponding to j            
            end
            M = size(obj.x,2);
            tMat = tMat(1:M*3,1:M*3);
            tVec = tVec(1:M*3);
            fprintf("GraphSLAM reduce completed\n")
        end
        
        function [mu] = solve(obj, tMat, tVec)
            % INPUT: tMat, tVec are info mat and vec with tilde
            % with tilde means only the poses, but not the map!
            % OUTPUT: mean value and coviance

            % GLOBAL VARIABLES
            numM = size(obj.m, 2);
            
            tic
            H_sparse = sparse(tMat);
            mu1 = H_sparse \ tVec;            
            toc
            %
%             cov = inv(tMat);
%             mu1 = cov * tVec;   % mu 0:t

            % compute mean value of map features
            mu2 = zeros(3 * numM,1);    % (with signature);

            for j = 1: numM
                tau_j = obj.tau{j};
                tau_j_idx = (3*(tau_j-1)+1):3*tau_j;   % tau(j): poses
                j_idx = find_iMat_idx(obj, j, 'm');  % j: map
                mu2(3*j-2:3*j) = obj.H(j_idx, j_idx) \ (obj.b(j_idx) + obj.H(j_idx, tau_j_idx) * tVec(tau_j_idx));
            end

            mu = [mu1; mu2];
            fprintf("GraphSLAM solve completed\n")
        end
        
        
        function plotgraph(obj)
            plot(obj.x(1,:),obj.x(2,:),'k');
            hold on
            plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);
            fprintf("Plot Pose Graph\n")
        end
        
        function run(obj)
            M = size(obj.x,2);
            N = size(obj.m,2);
            
            linearize(obj);
            [tMat, tVec] = reduce(obj);
            %[mu, cov] = solve(obj, tMat, tVec);
            [mu] = solve(obj, tMat, tVec);
            
            
            Mu = reshape(mu, 3, size(mu,1)/3);
            
            hold on
            plot(obj.x(1,:),obj.x(2,:),'b');    
            plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);
            plot(Mu(1,1:M)+obj.x(1,:),Mu(2,1:M)+obj.x(2,:),'r')
            plot(Mu(1,M+1:M+N)+obj.m(1,:),Mu(2,M+1:M+N)+obj.m(2,:),'r.','MarkerSize', 4)
%           plot(Mu(1,1:M),Mu(2,1:M),'r')
%           plot(Mu(1,M+1:M+N),Mu(2,M+1:M+N),'r.','MarkerSize', 4)
            

        end
     
        
        function idx = find_iMat_idx(obj, num, str)
            % get indices of the info matrix
            % num: can be a single num or a list
            % str: 'x': pose; 'm': maps

            % GLOBAL VARIABLES
            numX = size(obj.x, 2);

            if str == 'x'
                len = size(num(:), 1);
                if len == 1
                    idx = [3 * num - 2, 3 * num - 1, 3 * num];
                    return
                end

                list = reshape(num, [len, 1]);
                c1 = 3 * list - 2;
                c2 = c1 + 1;
                c3 = c2 + 1;
                c = [c1, c2, c3];
                c1 = c';
                idx = c1(:);

            end

            if str == 'm'
                bias = 3 * numX;
                len = size(num(:), 1);
                if len == 1
                    idx = [3 * num - 2, 3 * num - 1, 3 * num];
                    idx = idx + bias;
                    return
                end

                list = reshape(num, [len, 1]);
                c1 = 3 * list - 2;
                c2 = c1 + 1;
                c3 = c2 + 1;
                c = [c1, c2, c3];
                c1 = c';
                idx = c1(:);
                idx = idx + bias;
            end
        end
        
    end
end