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
       
       is_init = "no";
       
   end
   methods
       %% function area
       
       function r = EKF_init(obj, P_, Q_, R_, function_f_, function_jf_, function_h_, function_jh_)
           
           % matrix init
           obj.P = P_;
           obj.Q = Q_;
           obj.R = R_;
           
           % function init
           obj.function_f = function_f_;
           obj.function_h = function_h_;
           obj.function_jf = function_jf_;
           obj.function_jh = function_jh_;
           
           % init ok
           obj.is_init = "ok";
           r = obj.is_init;
       end
       
       function r = EKF_run(obj,x_pre_, u_, z_)
           if obj.is_init == "ok"
                x_size = size(x_pre_, 1);
                arguments = num2cell([x_pre_' u_']);
                F = obj.function_jf(arguments{:});
                H = obj.function_jh(arguments{:});
                
                state_hat_temp = obj.function_f(arguments{:});    % Prediction
                obj.P = F * obj.P * F' + obj.Q;
                K = obj.P * H' / (H*obj.P*H' + obj.R);
                Inno = z_ - obj.function_h(arguments{:});  % Innovation
                state_hat = state_hat_temp + K * Inno;   % Correction
                obj.P = (eye(x_size) - K*H) * obj.P * (eye(x_size) - K*H)' + K*obj.R*K';

                r = state_hat;
           else
               error("you must init class    : Call (filtering_init(obj, ...)");
           end
       end
   end
end







