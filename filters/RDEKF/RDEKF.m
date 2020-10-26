% paper:
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-07-28
% Extended Kalman Filter algorithm.



classdef RDEKF < handle
    properties
        
        % functions
        function_f;
        function_jf;
        function_h1;
        function_jh1;
        function_h2;
        function_jh2;
        function_h3;
        function_jh3;
        
        % covariance matrix
        P;
        % error matrix
        Q;
        R;
        % neighbor num
        nn;
        % data saving
        count = 1;
        x_appended;
        x_pre;
        % estimation error
        x_se;
        x_rmse;
        is_init = "no";
        first_run = 0;
    end
    methods
        %% function area
        function obj = RDEKF(P_, Q_, R_, function_f_, function_jf_, function_h1_, function_h2_, function_h3_, function_jh1_, function_jh2_, function_jh3_, init_, nn_)
            % matrix init
            obj.P = P_;
            obj.Q = Q_;
            obj.R = R_;
            
            % function init
            obj.function_f = function_f_;
            obj.function_jf = function_jf_;
            obj.function_h1 = function_h1_;
            obj.function_jh1 = function_jh1_;
            obj.function_h2 = function_h2_;
            obj.function_jh2 = function_jh2_;
            obj.function_h3 = function_h3_;
            obj.function_jh3 = function_jh3_;
            
            obj.nn = nn_;
            % x_hat appended
            obj.x_appended = zeros(size(init_,1), []);
            obj.x_appended(:,1) = init_(:);
            
            obj.x_pre = init_;
            
            % init ok
            obj.count = 1;
            obj.is_init = "ok";
        end
        
        function r = estimate3(obj, i_, u_, z_, adj_, pj_)
            if obj.is_init == "ok"
                x_size = size(obj.x_pre, 1);
                argument_f = num2cell([obj.x_pre' u_']);
                f_hat = obj.function_f(argument_f{:});    % Prediction
                F = obj.function_jf(argument_f{:});
                
                % make h_hat & H (jacobian)
                %                 argument_h = num2cell([f_hat' pj_']);
                %                 H = obj.function_jh(argument_h{:});
                find_neighbors = find(adj_(:,i_)==1);
                h_hat = zeros(2*obj.nn + 1,1);
                H = zeros(obj.nn*2+1, x_size);
                for i = 1:obj.nn
                    argument_h = num2cell([f_hat' pj_(1:2,find_neighbors(i))']);
                    h_hat(i,1) = obj.function_h1(argument_h{:});
                    H(i,:) = obj.function_jh1(argument_h{:});
                    h_hat(obj.nn+i,1) = obj.function_h2(argument_h{:});
                    H(obj.nn+i,:) = obj.function_jh2(argument_h{:});
                end
                argument_h = num2cell([f_hat' [0 0]]);
                h_hat(2*obj.nn + 1,1) = obj.function_h3(argument_h{:});
                H(2*obj.nn + 1,:) = obj.function_jh3(argument_h{:});
                
                obj.P = F * obj.P * F' + obj.Q;
                K = obj.P * H' / (H*obj.P*H' + obj.R);
                Inno = z_ - h_hat;  % Innovation
                state_hat = f_hat + K * Inno;   % Correction
                obj.P = (eye(x_size) - K*H) * obj.P * (eye(x_size) - K*H)' + K*obj.R*K';
                
                r = state_hat;
                obj.x_appended(:,obj.count) = r;
                obj.x_pre = r;
                obj.count = obj.count + 1;
            else
                error("you must init class    : Call (filtering_init(obj, ...)");
            end
        end
    end
end







