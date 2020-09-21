% paper: Distributed Multirobot Localization
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-09-07
% Extended Kalman filter, message exchange.
% Argumented system

classdef EKFDCL < handle
    properties
        Phi;            % function_f...(1)
        
        F;
        J_F;
        H;
        
        w;              % not used
        v;              % not used
        a;              % adjacency? or weight?
        P;              % covariance matrix
        Q;              % error covariance...(4)
        R;              % error matrix
        S;              % correlation term
        x_appended;     % data saving
        ns;             % neighbor size
        an;             % now robot num
        nx;             % state size
        nu;             % input size
        x_pre;          % previous state x; because of recursive algorithm.
        count = 1;      % The count variable are started from =2, beacause count==1 is for initialized state value.
        is_init = "no"; % initialized confirm
        % inner debuging
    end
    methods
        %% function area
        function obj = EKFDCL(P_, Q_, R_, F_, J_F_, a_, init_, agent_num_)
            % arguments
            % P_ initial cov
            % Q_ w cov
            % V_ y cov
            % R_ r cov
            % A B C system, measurement matrix
            % a_ adjacency graph (if i, j are setted 1 than, i obtain
            % measurement from j agent.
            % init_ initial state x
            % robot_num_ agent i robot
            if nargin == 0
                error("No params error. Initialized the class with (P,Q,R,A,B,C,D,init,algorithm_num)");
            else
                disp(nargin)
                % if k =1, initialize x_hat and initial P
                obj.P = P_;
                obj.Q = Q_;
                obj.R = R_;
                
                
                obj.F = F_;
                obj.J_F = J_F_;
                
                obj.a = a_; % adjacency
                obj.ns = size(a_,1); %neighbor size
                obj.an = agent_num_; % now robot num
                
                obj.x_appended = zeros(size(init_,1), []);
                obj.x_appended(:,1) = init_(:);
                obj.x_pre = init_;
                obj.count = 2;
                
                obj.nx = 3;
                obj.nu = 2;
                obj.is_init = "ok";
            end
            
        end
        
        function r = estimate_no_relative(obj, u_,y_)
           x_m = zeros(obj.nx * obj.an,1);
           a_F = zeros(obj.nx * obj.an, obj.nx * obj.an);
            for i = 1:obj.an
               start = 1+(i-1)*obj.nx;
               finish = obj.nx+(i-1)*obj.nx;
               startu = 1+(i-1)*obj.nu;
               finishu = obj.nu+(i-1)*obj.nu;
               arguments = num2cell([obj.x_pre(start:finish)' u_(startu:finishu)']);
%                disp(arguments)
               a_F(start:finish,start:finish) = obj.J_F(arguments{:});
               x_m(1+(i-1)*obj.nx:obj.nx+(i-1)*obj.nx) = obj.F(arguments{:});
            end
            obj.P = obj.P;
            obj.x_pre = x_m;
            r = x_m;
            obj.x_appended(:,obj.count) = r;
            obj.count = obj.count + 1;
        end
        
        function r = estimate(obj, agent_i, agent_j, u_, y_)
           % arguments
           % agent_i_index = i agent (obtain relative measurement)
           % agent_j_index = j agent target agent
           % relative_ = relative measurement from i to j (j-i)
           
           %prediction
           x_m = zeros(obj.nx * obj.an,1);
           a_F = zeros(obj.nx * obj.an, obj.nx * obj.an);
           
           
           
           for i = 1:obj.an
               start = 1+(i-1)*obj.nx;
               finish = obj.nx+(i-1)*obj.nx;
               startu = 1+(i-1)*obj.nu;
               finishu = obj.nu+(i-1)*obj.nu;
               arguments = num2cell([obj.x_pre(start:finish)' u_(startu:finishu)']);
%                disp(arguments)
               a_F(start:finish,start:finish) = obj.J_F(arguments{:});
               x_m(1+(i-1)*obj.nx:obj.nx+(i-1)*obj.nx) = obj.F(arguments{:});
           end
           
           obj.P = a_F * obj.P * a_F' + obj.Q;
           
           I3 = eye(3);
           obj.H = zeros(obj.nx, obj.nx * obj.an);
           obj.H(1:3,1+(agent_i-1)*obj.nx:(agent_i-1)*obj.nx+obj.nx) = -I3;
           obj.H(1:3,1+(agent_j-1)*obj.nx:(agent_j-1)*obj.nx+obj.nx) = I3;
           
           obj.S = obj.H*obj.P*obj.H'+obj.R;
           K = obj.P * obj.H' / obj.S;
           
           predict = obj.H * x_m;
           
           inno = y_ - predict;
           
           r = x_m + K * inno;
           
           obj.P = (eye(obj.nx*obj.an) - K * obj.H) * obj.P * (eye(obj.nx*obj.an) - K*obj.H)' + K * obj.R*K';
           
           %update
           
            obj.x_pre = r;
            obj.x_appended(:,obj.count) = r;
            obj.count = obj.count + 1;
        end
        
    end
end



