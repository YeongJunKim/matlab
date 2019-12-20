classdef FF_EKF_UFIR < handle
   properties 
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       
       % x_pre_
       x_pre = [];
       % measurement
       z = [];
       % control input
       u = [];
       % arguments
       arguments_f = [];
       arguments_h = [];
       % esitimation x
       x_hat_KF = [];
       x_hat_UFIR = [];
       x_hat = [];
       % likelihood
       likelihood = zeros(2,1);
       % likelihood weight
       weight = zeros(2,1);
       % Markov transition probability matrix
       markov = zeros(2,2);
       
       
       % P Q R
       P;
       Q;
       R;
       
       % horizon size
       h_size;
       % x size
       x_size;
       % u size
       u_size;
       % z size
       z_size;
        
       % system matrix set
       f_tilde;
       % measurement update matrix set
       h_tilde;
       % measurement set
       y_tilde;
       % control input set
       u_tilde;
       % state set
       x_tilde;
       % noise set
       w_tilde;
       v_tilde;
       
       
       % init state
       init_state;
       
       
       
       % counting
       count = 0;
       
       is_init = "no";
       
       % data saving
       
       x_appended;
       
       x_appended_EKF;
       x_appended_UFIR;
       
   end
   methods
       %% function area
       function r = FF_EKF_UFIR_init(obj, h_size_, x_size_, z_size_, u_size_, function_f_, function_jf_, function_h_, function_jh_, init_state_, P_, Q_, R_, weight_, markov_)
           % size mapping
           obj.h_size = h_size_;
           obj.x_size = x_size_;
           obj.u_size = u_size_;
           obj.z_size = z_size_;

           % array init
           obj.f_tilde = zeros(x_size_, x_size_, h_size_);
           obj.h_tilde = zeros(z_size_, x_size_, h_size_);
           obj.y_tilde = zeros(z_size_, h_size_);
           obj.u_tilde = zeros(x_size_, h_size_);
           
           % function init
           obj.function_f = function_f_;
           obj.function_h = function_h_;
           obj.function_jf = function_jf_;
           obj.function_jh = function_jh_;
           
           % P Q R
           obj.P = P_;
           obj.Q = Q_;
           obj.R = R_;
           
           % weight
           obj.weight = weight_;
           obj.markov = markov_;
           
           % count init
           obj.count = 1;
           
           % state init
           obj.init_state = init_state_;
           
           % data saving
           obj.x_appended = zeros(x_size_, 1000);
           
           % init ok
           obj.is_init = "ok";
           r = obj.is_init;
       end
       
       function r = FF_EKF_UFIR_run(obj, x_pre_, u_, z_)
           if obj.is_init == "ok"
               obj.arguments_f = num2cell([x_pre_' u_']);
               f_hat = obj.function_f(obj.arguments_f{:});
               obj.arguments_h = num2cell(f_hat');
               h_hat = obj.function_h(obj.arguments_h{:});
               
               F =  obj.function_jf(obj.arguments_f{:});
               H =  obj.function_jh(obj.arguments_h{:});
               
               % accumulating array
               obj.f_tilde(:,:,1:obj.h_size - 1) = obj.f_tilde(:,:,2:obj.h_size);
               obj.f_tilde(:,:,obj.h_size) = F;
               obj.h_tilde(:,:,1:obj.h_size-1) = obj.h_tilde(:,:,2:obj.h_size);
               obj.h_tilde(:,:,obj.h_size) = H;
               obj.y_tilde(:,1:obj.h_size-1) = obj.y_tilde(:,2:obj.h_size);
               obj.y_tilde(:,obj.h_size) = z_ - (h_hat - H * f_hat);
               obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
               obj.u_tilde(:,obj.h_size) = f_hat - F * x_pre_;
               obj.x_tilde(:,1:obj.h_size-1) = obj.x_tilde(:,2:obj.h_size);
               obj.x_tilde(:,obj.h_size) = F * x_pre_ + obj.w_tilde;
               
               obj.x_pre = x_pre_;
               obj.u = u_;
               obj.z = z_;
          
                              
               if obj.count > obj.h_size
                     r1 = FF_UFIR_main(obj);
                     r2 = FF_EKF_main(obj);
                     
               else
                     r = FF_EKF_main(obj);
               end
               
               obj.count = obj.count + 1;
               
           else
               error("you must init class    : Call (@@_init(obj, ...)");
           end
       end
   end
end



function r = FF_UFIR_main(obj)
    x_size = size(x_pre_, 1);
    state_hat_temp = obj.function_f(obj.arguments_f{:});    % Prediction
    obj.P = F * obj.P * F' + obj.Q;
    K = obj.P * H' / (H*obj.P*H' + obj.R);
    Inno = z_ - obj.function_h(obj.arguments_h{:});  % Innovation
    state_hat = state_hat_temp + K * Inno;   % Correction
    obj.P = (eye(x_size) - K*H) * obj.P * (eye(x_size) - K*H)' + K*obj.R*K';

    r = state_hat;
    
    if obj.first_run == 0
        obj.first_run = 1;
        obj.x_appended_EKF = zeros(x_size, 1000);
    end
    obj.x_appended_EKF(:,obj.count) = r;
end

function r = FF_EKF_main(obj)

end


















