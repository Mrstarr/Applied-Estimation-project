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
        N_ldm   % Number of landmarks
        N_X     % Number of states 
        N_Z     % Number of measurement
    end  % properties set private
    
    properties (Dependent = true)
        
        % pose    % Poses of all nodes
    end  % properties dependent
    
    methods

        % initialize the pose graph
        function initialize(obj,time,TLsr,u,zt)
            %I = round(size(u,2)/20);
            %I = 10;
            global I 
            ILsr = size(zt,2);
            
            % u & x
            obj.u = zeros(3, I-1);
            obj.x = zeros(3,I);
            % get control data
            obj.u(1,:) = u(1,1:I-1); % motion speed
            obj.u(2,:) = u(2,1:I-1); % steering angle
            obj.u(3,:) = diff(time(1,1:I)); % time difference           
            
            global x0
            obj.x(:,1) = [0;0;0];
            
            % initialization
            j = 1;
            obj.N_ldm = 0;
            obj.N_Z = 0;
            obj.N_X = 1;  % x0 already exists
            obj.z.m = []; % measurement 
            obj.z.idx = []; % time index i
            obj.z.c = []; % correspondence
            
            % input all control info
            for i = 1:size(obj.u,2)
                obj.x(:,i+1) = predict(obj.x(:,i),obj.u(:,i));
                obj.N_X = obj.N_X + 1; 
                % when new landmarks detected
                if j <= ILsr
                    if TLsr(j) < time(i)
                        % add new measurement
                        n_zt = size(zt{j},2);
                        obj.z.m = [obj.z.m zt{j}];  
                        obj.z.idx = [obj.z.idx i * ones(1,n_zt)]; 
                        obj.z.c = [obj.z.c (obj.N_ldm + 1):obj.N_ldm + n_zt]; 
                        obj.N_Z = obj.N_Z + n_zt;
                        
                        % add new landmark
                        ldm = add_landmark(obj.x(:,i),zt{j});   
                        obj.m = [obj.m ldm];
                        for q = (obj.N_ldm + 1) : (obj.N_ldm + size(zt{j},2))
                            obj.tau{q} = i; % tau{q}: list of index of x where detects landmark q 
                        end
                        obj.N_ldm = obj.N_ldm + n_zt;
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
            R = noise.R;    % motion error
            Q = noise.Q;    % measurement error
            % initialize H and b
            obj.n_node = obj.N_X + obj.N_ldm;
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
            for j = 1:obj.N_Z
                
                c = obj.z.c(j);
                t = obj.z.idx(j);
                dlt = [obj.m(1,c) - obj.x(1,t);
                         obj.m(2,c) - obj.x(2,t)];
                q = dlt'*dlt;
                z_hat = [sqrt(q);atan2(dlt(2),dlt(1))-obj.x(3,t)+pi/2;obj.m(3,c)];
                H_1 = 1/q * [ -sqrt(q) * dlt(1) -sqrt(q) * dlt(2) 0;
                              dlt(2)  -dlt(1)  -q;
                              0 0 0];
                H_2 = 1/q * [ sqrt(q) * dlt(1) sqrt(q) * dlt(2) 0;
                              -dlt(2)  dlt(1)  0;
                              0 0 q];

                % update information matrix       
                t_idx = (3*(t-1)+1):3*t;
                j_idx = (3*(c-1)+1+3*obj.N_X):(3*c+3*obj.N_X);
                H_11 = H_1' / Q * H_1;
                H_12 = H_1' / Q * H_2;
                H_21 = H_2' / Q * H_1;
                H_22 = H_2' / Q * H_2;
                b1 = H_1' / Q * (obj.z.m(:,j) - z_hat + [H_1 H_2]*[obj.x(:,t); obj.m(:,c)]);
                b2 = H_2' / Q * (obj.z.m(:,j) - z_hat + [H_1 H_2]*[obj.x(:,t); obj.m(:,c)]);

                obj.H(t_idx,t_idx) = obj.H(t_idx,t_idx) + H_11;
                obj.H(t_idx,j_idx) = obj.H(t_idx,j_idx) + H_12;
                obj.H(j_idx,t_idx) = obj.H(j_idx,t_idx) + H_21;
                obj.H(j_idx,j_idx) = obj.H(j_idx,j_idx) + H_22;
                obj.b(t_idx) = obj.b(t_idx) + b1;
                obj.b(j_idx) = obj.b(j_idx) + b2;
            end
            
            fprintf("Linearization Complete\n")
            fprintf("Size of Information matrix: %d * %d\n", size(obj.H))
        end       
        
        function [tMat, tVec] = reduce(obj)
            tMat = obj.H;
            tVec = obj.b;
            
            for j =1:obj.N_ldm
                for g = obj.tau{j}
                    % find corresponding index of j and x
                    idx_x = find_iMat_idx(obj,g,"x");
                    idx_j = find_iMat_idx(obj,j,"m");
                    % reduce H and b
                    tVec(idx_x) = tVec(idx_x)- tMat(idx_x,idx_j)/tMat(idx_j,idx_j)*tVec(idx_j);
                    tMat(idx_x,idx_x) = tMat(idx_x,idx_x)- tMat(idx_x,idx_j)/tMat(idx_j,idx_j)*tMat(idx_j,idx_x);
                end
                % remove rows/columns corresponding to j            
            end
            tMat = tMat(1:obj.N_X*3,1:obj.N_X*3);
            tVec = tVec(1:obj.N_X*3);
            fprintf("GraphSLAM reduce completed\n")
        end
        
         function [mu, cov] = solve(obj, tMat, tVec)
            % INPUT: tMat, tVec are info mat and vec with tilde
            % with tilde means only the poses, but not the map!
            % OUTPUT: mean value and coviance

            % GLOBAL VARIABLES
            numM = obj.N_ldm;
            
            tic
            H_sparse = sparse(tMat);
            cov = inv(H_sparse);
            mu1 = H_sparse \ tVec;            
            toc
            %
%             cov = inv(tMat);
%             mu1 = cov * tVec;   % mu 0:t

            % compute mean value of map features
            mu2 = zeros(3 * numM,1);    % (with signature);

            for j = 1: numM
                tau_j = obj.tau{j};
                tau_j_idx = find_iMat_idx(obj, tau_j,'x');   % tau(j): poses
                j_idx = find_iMat_idx(obj, j, 'm');  % j: map
                mu2(3*j-2:3*j) = obj.H(j_idx, j_idx) \ (obj.b(j_idx) + obj.H(j_idx, tau_j_idx) * tVec(tau_j_idx));
            end

            mu = [mu1; mu2];
            fprintf("GraphSLAM solve completed\n")
        end
        
        function prob = correspondence_test (obj, mu, covar, j, k)
            % TARGET: get mean vector mu and the path covariance covar(0~t)
            % from GraphSLAM_solve
            % INPUT: iVec not needed. deleted.
            % OUTPUT: the posterior prob of mj & mk are the same

            jk = find_iMat_idx(obj, [j, k], 'm');          % mat idx: feature j&k 
            tau_jk_num = union(obj.tau{j},obj.tau{k});                    % the tau poses of j
            tau_jk = find_iMat_idx(obj, tau_jk_num, 'x');  % mat idx: tau-poses

            % tau_set() should be implemented in the graph class
            % the set of poses ? (j, k) at which the robot observed feature j and k

            iMat_jk = obj.H(jk, jk);
            iMat_jkt = obj.H(jk, tau_jk);
            iMat_tjk = obj.H(tau_jk, jk);
            cov_jk = covar(tau_jk, tau_jk);

            % marginalized info
            iMat_margin = iMat_jk - iMat_jkt * cov_jk * iMat_tjk;
            iVec_margin = iMat_margin * mu(jk) ;

            % info of the difference variable
            iMat_diff = [eye(3), -eye(3)] * iMat_margin * [eye(3); -eye(3)];
            iVec_diff = [eye(3), -eye(3)] * iVec_margin;
            %mu_diff = iMat_diff \ iVec_diff;
            mu_diff = [eye(3), -eye(3)] * mu(jk);

            %prob = sqrt(det(iMat_diff))/sqrt(2*pi) * ...
            %    exp(-0.5 * mu_diff' / iMat_diff * mu_diff);
            prob = 1/sqrt(2*pi) * ...
                exp(-0.5 * abs(mu_diff' / iMat_diff * mu_diff));

        end
        
        function plotgraph(obj)
            plot(obj.x(1,:),obj.x(2,:),'k');
            hold on
            plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);
            fprintf("Plot Pose Graph\n")
        end
        
        function batch_association(obj,mu,cov)
            N = obj.N_ldm;
            tic
            for m_j = 1:N
                for m_k = m_j+1:N
                    if obj.z.c(m_k) < m_k
                        continue
                    end
                    
                    P_jk = correspondence_test(obj,mu,cov,m_j,m_k);
                    %fprintf("P between %d and %d is %d\n",m_j,m_k,P_jk);
                    if P_jk > 0.35
                        obj.z.c(m_k) = m_j;
                    end
                end
                if mod(m_j,100) == 0
                    fprintf("Batch Association progress: %d\n",m_j)
                end
            end
            toc         
        end
        
        function update_map(obj)
            N = size(obj.z.m,2);
            idx = 1;    % new landmark index 
            m_temp =[]; % new map
            
            % create look up table
            for i =1:N
                if obj.z.c(i) == i 
                    idx_lookup(i) = idx; % shrink map by look up table
                    idx = idx + 1; 
                    m_temp = [m_temp obj.m(:,i)]; % update new map
                else 
                    idx_lookup(i) = idx_lookup(obj.z.c(i));
                end
            end
            obj.m = m_temp;
            obj.N_ldm = size(obj.m,2);
            
            % shrink map(index related to map)
            % 1. shrink correspondence
            
            for i = 1:obj.N_Z
                c_idx = obj.z.c(i);
                obj.z.c(i) = idx_lookup(c_idx);
            end
            
            % 2. shrink tau
            tau_temp = cell(1,obj.N_ldm);
            for i = 1:size(obj.tau,2)
                tau_idx = idx_lookup(i);
                tau_temp{tau_idx} = union(tau_temp{tau_idx},obj.tau{i});
            end
            obj.tau = tau_temp;
            
            fprintf("Map Update Complete\n")        
        end           
            
        function run(obj, verbose)
            close all
            linearize(obj);
            [tMat, tVec] = reduce(obj);
            [mu, cov] = solve(obj, tMat, tVec);
            if verbose == 2
                batch_association(obj,mu,cov)
                update_map(obj)
                hold on
                plot(obj.x(1,:),obj.x(2,:),'k');    
                plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);
                close all 
                
                linearize(obj);
                [tMat, tVec] = reduce(obj);
                [mu, cov] = solve(obj, tMat, tVec);
                 Mu = reshape(mu, 3, size(mu,1)/3);
                 hold on
                 plot(Mu(1,1:obj.N_X),Mu(2,1:obj.N_X),'k');
                 plot(Mu(1,obj.N_ldm + 1:obj.N_ldm +obj.N_X),Mu(2,obj.N_ldm +1:obj.N_ldm +obj.N_X),'MarkerSize', 4);
                 close all
            end
            
            if verbose == 3
                j = 1;
                while j < 8
                    batch_association(obj,mu,cov)
                    update_map(obj)
                    
                    hold on
                    plot(obj.x(1,:),obj.x(2,:),'k');    
                    plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);
                    close all 
                    linearize(obj);
                    [tMat, tVec] = reduce(obj);
                    [mu, cov] = solve(obj, tMat, tVec);
                    
                    j = j + 1;
                end
            end
            Mu = reshape(mu, 3, size(mu,1)/3);
            M = size(obj.x,2);
            N = size(obj.m,2);
            obj.x =  Mu(:,1:M);
            obj.m =  Mu(:,M+1:M+N);


%             Mu = reshape(mu, 3, size(mu,1)/3);
            state = obj.x;
            landmark = obj.m;
            save("state_1.mat",'state')
            save("landmark_1.mat",'landmark')

            
            plot(obj.x(1,:),obj.x(2,:),'k');
            hold on
            plot(obj.m(1,:),obj.m(2,:),'b.','MarkerSize', 4);

%             plot(Mu(1,1:M)+obj.x(1,:),Mu(2,1:M)+obj.x(2,:),'r')
%             plot(Mu(1,M+1:M+N)+obj.m(1,:),Mu(2,M+1:M+N)+obj.m(2,:),'r.','MarkerSize', 4)
%             plot(Mu(1,1:M),Mu(2,1:M),'r')
%             plot(Mu(1,M+1:M+N),Mu(2,M+1:M+N),'r.','MarkerSize', 4)
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