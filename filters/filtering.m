classdef filtering < handle
   properties 
       %% variable area
       
       filter_type;
       
       f,h,jac_f,jac_h;
       x_hat = [];
       x_hat_pre = [];
       
       z = [];
       u = [];
       
       alpha = 0;
       h_size = 0;
       isset = 0;
       time_step = 1;
       
        % function [x_hat] = PEFFME_main(f, h, jac_f, jac_h, x_pre, z, u, alpha, h_size)
        % function [x_hat] = EKF_main(f, h, jac_f, jac_h, x_pre, z, u, P, Q, R) 
   end
   methods
       %% function area
       function r = filtering_init(obj,method, f,h,jac_f,jac_h)
           addpath(genpath());
           obj.f = f;
           obj.h = h;
           obj.jac_f = jac_f;
           obj.fac_h = jac_h;
           if obj.isset == 0
               if method == EKF
                    obj.filter_type = @EKF_main;
                    r = "OK";
               elseif method == FIR
                    obj.filter_type = @FIR_main;
                    r = "OK";
               elseif method == KF
                    obj.filter_type = @KF_main;
               r = "OK";
               elseif method == PEFFME
                    obj.filter_type = @PEFFME_main;
               r = "OK";
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
       function x_hat = filtering_run(obj, measurement)
           if obj.isset == 1
           else
               disp("NOT SETED FILTER");
           end
       end
       
   end
end