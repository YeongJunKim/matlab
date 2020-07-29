% paper: Distributed Kalman filter for cooperative localization with integrated measurements
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-07-28
% There are two algorithm, one is the algorithm1 and other is algorithm2.
% Each algorithm weighted more in absoluted measurement or relative measurement.
% So I will implement these two algorithm in one class function.
% Moreover this class is fitted in linear system.
% Utilezed code: https://www.github.com/YeongJunKim/localization_sim

classdef DKFCL < handle
    properties
        algorithm_num   % algorithm type: can be define 1, 2, or 3. 3 is the general kalman filter L to 0.
        A;              % A matrix ...(1)
        B;              % B matrix ...(1)
        C;              % absolute measurement matrix Cx + v ...(2)
        D;              % relative measurement matrix D(x_diff) + r ...(3)
        w;              % not used
        v;              % not used
        a;              % adjacency? or weight?
        P;              % covariance matrix
        Q;              % error matrix
        V;              % error matrix
        R;              % error matrix
        x_appended;     % data saving
        ns;             % neighbor size
        rn;             % now robot num
        nx;             % state size
        x_pre;          % previous state x; because of recursive algorithm.
        count = 1;      % The count variable are started from =2, beacause count==1 is for initialized state value.
        is_init = "no"; % initialized confirm
        % inner debuging
        y               % debuging y
        z               % debuging z
        x_neigh         % debuging x_neighbor's
    end
    methods
        %% function area
        function obj = DKFCL(P_, Q_, V_, R_, A_, B_, C_, D_,  a_, init_, algorithm_num_, robot_num_)
            % arguments
            % P_ initial cov
            % Q_ w cov
            % V_ y cov
            % R_ r cov
            % A B C D system, measurement matrix
            % a_ adjacency weighted graph
            % init_ initial state x
            % algorithm_num_ can be binary 0 or 1
            % robot_num_ agent i robot
            if nargin == 0
                error("No params error. Initialized the class with (P,Q,R,A,B,C,D,init,algorithm_num)");
            else
                
                
                % Algorithm 1 and 2 are saved in class but $FILTER$_run function return specific algorithm's output.
                % Choice which algorithm do you want.
                obj.algorithm_num = algorithm_num_;
                
                % if k =1, initialize x_hat and initial P
                obj.P = P_;
                obj.Q = Q_;
                obj.V = V_;
                obj.R = R_;
                obj.A = A_;
                obj.B = B_;
                obj.C = C_;
                obj.D = D_;
                obj.a = a_; % adjacency
                obj.ns = size(a_,1); %neighbor size
                obj.rn = robot_num_; % now robot num
                
                
                obj.x_appended = zeros(size(init_,1), []);
                obj.x_appended(:,1) = init_(:);
                obj.x_pre = init_;
                obj.count = 2;
                
                obj.nx = size(init_,1);
                
                obj.is_init = "ok";
                
                if(obj.algorithm_num ~= 1 && obj.algorithm_num ~= 2)
                    error("Algorithm number can be the number as 1 or 2.");
                end
                
            end
        end
        
        function r = DKFCL_obtain_x_prior(obj)
            r = obj.A * obj.x_pre;
        end
        
        function r = DKFCL_obtain_y_prior(obj)
            r = obj.A * obj.P * obj.A' + obj.B * obj.Q * obj.B';
        end
        
        function r = DKFCL_run(obj, u_, y_, z_, x_neighbor, p_neighbor)
            % arguments
            % u_ control input
            % y_ absolute measurement
            % z_ relative measurement
            % x_neighbor neighbor's prior_x size == (nx, ns)
            % p_neighbor neighbor's prior_p size == (nx, nx, ns)
            obj.y = y_;
            obj.z = z_;
            obj.x_neigh = x_neighbor;
            s = obj.is_init;
            if strcmp(s,"ok") == 0
                error("you must init class DKFCL(arg..)");
            end
            obj.x_pre(3:4) = u_(1:2);
            if(obj.algorithm_num == 1)
                % 2. sample integrated measurements y_{i,k} and z_{i,j,k}
                % 3. predict target estimates
                prior_x = obj.A * obj.x_pre;
                prior_p = obj.A * obj.P * obj.A' + obj.B * obj.Q * obj.B';
                % 4. send bar{x}_{i,k} to target j,
                % 5. receive bar{x}_{j,k} and bar{P}_{j,k} from target j,
                % 6. calculate gain matrices
                % algorithm 1
                temp = zeros(2,2);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(obj.rn,j) ~= 0)
                        temp = temp +  obj.a(i,j)^2 * obj.R;
                    end
                end
                K = prior_p * obj.C' * (obj.C * prior_p * obj.C' + obj.V)^-1;
                L = (prior_p * obj.D' - K * obj.C * prior_p * obj.D') * (obj.D * prior_p * obj.D' + temp)^-1;
                
                % 7. update target estimates
                z_hat = zeros(2, obj.ns);
                for j = 1:obj.ns
                    if(obj.a(obj.rn,j) ~= 0)
                        z_hat(:,j) = obj.D * (prior_x - x_neighbor(:,j));
                    end
                end
                temp = zeros(2,1);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(i,j) ~= 0)
                        temp = temp + obj.a(i,j) * (z_(1:2,j) - z_hat(:,j));
                    end
                end
                obj.x_appended(:,obj.count) = prior_x + K * (y_ - obj.C * prior_x) + L * temp;
                obj.x_pre(:) = obj.x_appended(:,obj.count);
                size_p = size(obj.P,1);
                temp1 = zeros(2,4);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(obj.rn,j) ~= 0)
                        temp1 = temp1 + (obj.a(i,j)^2 * obj.R * L');
                    end
                end
                temp2 = zeros(obj.nx,obj.nx);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(obj.rn,j) ~= 0)
                        temp2 = temp2 + obj.a(i,j)^2 * p_neighbor(:,:,j)* obj.D' * L';
                    end
                end
                obj.P = (eye(size_p) - K * obj.C - L * obj.D) * prior_p * (eye(size_p) - K * obj.C - L * obj.D)'+ K * obj.V * K' + L * temp1 + L * obj.D * temp2;
            elseif(obj.algorithm_num == 2)
                % 2. sample integrated measurements y_{i,k} and z_{i,j,k}
                % 3. predict target estimates
                prior_x = obj.A * obj.x_pre;
                prior_p = obj.A * obj.P * obj.A' + obj.B * obj.Q * obj.B';
                % 4. send bar{x}_{i,k} to target j,
                % 5. receive bar{x}_{j,k} and bar{P}_{j,k} from target j,
                % 6. calculate gain matrices
                % algorithm 2
                temp = zeros(2,2);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(i,j) ~= 0)
                        temp = temp +  obj.a(i,j)^2 * obj.R;
                    end
                end
                L = prior_p * obj.D' * (obj.D * prior_p * obj.D' + temp)^-1;
                K = (prior_p * obj.C' - L * obj.D * prior_p * obj.C') * (obj.C * prior_p * obj.C' + obj.V)^-1;
                % 7. update target estimates
                z_hat = zeros(2, obj.ns);
                for j = 1:obj.ns
                    if(obj.a(obj.rn,j) ~= 0)
                        z_hat(:,j) = obj.D * (prior_x - x_neighbor(:,j));
                    end
                end
                temp = zeros(2,1);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(i,j) ~= 0)
                        temp = temp + obj.a(i,j) * (z_(1:2,j) - z_hat(:,j));
                    end
                end
                
                obj.x_appended(:,obj.count) = prior_x + K * (y_ - obj.C * prior_x) + L * temp;
                obj.x_pre(:) = obj.x_appended(:,obj.count);
                size_p = size(obj.P,1);
                temp1 = zeros(2,4);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(i,j) ~= 0)
                        temp1 = temp1 + obj.a(i,j)^2 * obj.R * L';
                    end
                end
                temp2 = zeros(4,4);
                for j = 1:obj.ns
                    i = obj.rn;
                    if(obj.a(i,j) ~= 0)
                        temp2 = temp2 + obj.a(i,j)^2 * p_neighbor(:,:,j)* obj.D' * L';
                    end
                end
                obj.P = (eye(size_p) - K * obj.C - L * obj.D) * prior_p * (eye(size_p) - K * obj.C - L * obj.D)'+ K * obj.V * K' + L * temp1 + L * obj.D * temp2;
                
            elseif(ojb.algorithm_num == 3)
                % 2. sample integrated measurements y_{i,k} and z_{i,j,k}
                % 3. predict target estimates
                prior_x = obj.A * obj.x_pre;
                prior_p = obj.A * obj.P * obj.A' + obj.B * obj.Q * obj.B';
                % 4. send bar{x}_{i,k} to target j,
                % 5. receive bar{x}_{j,k} and bar{P}_{j,k} from target j,
                % 6. calculate gain matrices
                % algorithm 1
                temp = zeros(2,2);
                for j = 1:obj.ns
                    i = obj.rn;
                    temp = temp +  obj.a(i,j)^2 * obj.R;
                end
                K = prior_p * obj.C' * (obj.C * prior_p * obj.C' + obj.V)^-1;
                % 7. update target estimates
                z_hat = zeros(2, obj.ns);
                for j = 1:obj.ns
                    if(obj.a(obj.rn,j) ~= 0)
                        z_hat(:,j) = obj.D * (prior_x - x_neighbor(:,j));
                    end
                end
                temp = zeros(2,1);
                for j = 1:obj.ns
                    i = obj.rn;
                    temp = temp + obj.a(i,j) * (z_(1:2,j) - z_hat(:,j));
                end
                obj.x_appended(:,obj.count) = prior_x + K * (y_ - obj.C * prior_x)  * temp;
                obj.x_pre(:) = obj.x_appended(:,obj.count);
            end
            
            r = obj.x_appended(:,obj.count);
            obj.count = obj.count + 1;
        end
        
    end
end







