% paper:
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-11-13
% Unscented Kalman Filter algorithm.


classdef UKF
    %UKF 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    properties
        % for simple UKF
        function_f;
        function_jf;
        function_h;
        function_jh;
        % add more for relative_measurement
        function_h2;
        function_jh2;
        function_h3;
        function_jh3;
        P; Q; R;
        count = 1;
        x_appended;
        x_pre;
        x_remse;
        is_init = "no";
        first_run = 0;
    end
    
    methods
        function obj = UKF(P_, Q_, R_, f_, jf_, h_, jh_, init_)
            % P Q R
            obj.P = P_; obj.Q = Q_; obj.R = R_;
            % function
            obj.function_f = f_;
            obj.function_h = h_;
            obj.function_jf = jf_;
            obj.function_jh = jh_;
            % init
            obj.x_appended = zeros(size(init_,1), []);
            obj.x_appended(:,1) = init_(:);
            obj.x_pre = init_;
            obj.count = 1;
            obj.is_init = "ok";
            r = obj.is_init;
        end
        function r = estimate(obj, u_, z_)
            
        end
        function r = estimate_relative(obj, i_, u_, z_, adj_, pj_)
            x_size = size(obj.x_pre, 1);
            arguments_f = num2cell([obj.x_pre' u']);
            f_hat = obj.function_f(argument_f{:});
            F = obj.function_jh(argument_f{:});
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 이 메서드의 요약 설명 위치
            %   자세한 설명 위치
            outputArg = obj.Property1 + inputArg;
        end
    end
end

