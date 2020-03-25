classdef EKF < handle
   properties 
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       
       % covariance matrix
       P;
       % error matrix
       Q;
       R;
       
       % data saving
       count = 1;
       x_appended;
       x_pre;
       
       is_init = "no";
       first_run = 0;
   end
   methods
       %% function area
       
       function r = EKF_init(obj, P_, Q_, R_, function_f_, function_jf_, function_h_, function_jh_, init_)
           
           % matrix init
           obj.P = P_;
           obj.Q = Q_;
           obj.R = R_;
           
           % function init
           obj.function_f = function_f_;
           obj.function_h = function_h_;
           obj.function_jf = function_jf_;
           obj.function_jh = function_jh_;
           
           % x_hat appended
           obj.x_appended = zeros(size(init_,1), 1000);
           obj.x_appended(:,1) = init_(:);
           
          obj.x_pre = init_;
           
           % init ok
           obj.count = 1;
           obj.is_init = "ok";
           r = obj.is_init;
       end
       
       function r = EKF_run(obj,x_pre_, u_, z_)
           if obj.is_init == "ok"
                x_size = size(obj.x_pre, 1);
                arguments = num2cell([obj.x_pre' u_']);
                F = obj.function_jf(arguments{:});
                H = obj.function_jh(arguments{:});
                
                state_hat_temp = obj.function_f(arguments{:});    % Prediction
                obj.P = F * obj.P * F' + obj.Q;
                K = obj.P * H' / (H*obj.P*H' + obj.R);
                Inno = z_ - obj.function_h(arguments{:});  % Innovation
                state_hat = state_hat_temp + K * Inno;   % Correction
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







