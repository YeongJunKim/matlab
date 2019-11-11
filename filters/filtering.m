classdef filtering < handle
   properties 
       %% variable area
       method;
       filter_type;
       
       f,h,jac_f,jac_h;
       A, H;
       P, Q, R;
       x_hat = [];
       x_hat_pre = [];
       
       z = [];
       u = [];
       
       h_size = 0;
       isset = 0;
       time_step = 1;
       
   end
   methods
       %% function area
       function r = filtering_init(obj,method_, f,h,jac_f,jac_h)
           
           obj.method = method_;
           obj.f = f;
           obj.h = h;
           obj.jac_f = jac_f;
           obj.jac_h = jac_h;
           if obj.isset == 0
               if obj.method == "EKF"
                    obj.filter_type = @EKF_main;
                    r = "OK";
               elseif obj.method == "FIR"
                    obj.filter_type = @FIR_main;
                    r = "OK";
               elseif obj.method == "KF"
                    obj.filter_type = @KF_main;
                    r = "OK";
               elseif obj.method == "PEFFME"
                    obj.filter_type = @PEFFME_main;
                    r = "OK";
               %% if there is no filter you want, you can add filter under lines.
               %  elseif method == your filter
               else
                    r = "ERROR";
               end
               
               if r == "OK"
                   obj.isset = 1;
               else
                   obj.isset = 0;
               end
           end
       end
       function r = filtering_set_h_size(obj, h_size_)
           obj.h_size = h_size_;
           r = "OK";
       end
       function r = filtering_set_PQR(obj, P_, Q_, R_)
          obj.P = P_;
          obj.Q = Q_;
          obj.R = R_;
          r = "OK";
       end
       function x_hat = filtering_run(obj, u_, z_, alpha_)
           if obj.method == "EKF"              
               x_hat = obj.filter_type(obj.f, obj.h, obj.jac_f, obj.jac_h, obj.x_hat_pre, z_, u_, obj.P, obj.Q, obj.R);
               obj.x_hat_pre = x_hat;
               obj.time_step = obj.time_step + 1;
           elseif obj.method == "FIR"
               
               obj.time_step = obj.time_step + 1;
           elseif obj.method == "KF"
                    xp = A*x;
                    Pp = A*P*A'+Q;
                    K = Pp*H'*inv(H*Pp*H'+R);
                    x = xp+K*(z-H*xp);
                    P = Pp-K*H*Pp;  
               obj.time_step = obj.time_step + 1;
           elseif obj.method == "PEFFME"
               x_hat = obj.filter_type(obj.f, obj.h, obj.jac_f, obj.jac_h, obj.x_hat_pre, z_, u_, alpha_, obj.h_size);
               obj.x_hat_pre = x_hat;
               obj.time_step = obj.time_step + 1;
           else
               disp("NOT SETED FILTER");
           end
       end
       
   end
end