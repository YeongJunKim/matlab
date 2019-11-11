classdef FIR_colson < handle
   properties 
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       
       % measurement
       z = [];
       % control input
       u = [];
       % esitimation x
       x_hat = [];
       
       % horizon size
       h_size;
       % x size
       x_size;
       % u size
       u_size;
       % z size
       z_size;
        
       % system matrix set
       F_array;
       % measurement update matrix set
       H_array;
       % measurement set
       y_tilde;
       % control input set
       u_tilde;
       
       % z_
       z_;
       % u_
       u_;
       
       % counting
       count = 0;
       
       is_init;
       
   end
   methods
       %% function area
       
       function r = FIR_colson_init(obj, h_size_, x_size_, z_size_, u_size_, function_f_, function_jf_, function_h_, function_jh_)
%        function r = FIR_colson_init(obj, h_size_, x_size_, z_size_, u_size_)
           % size mapping
           obj.h_size = h_size_;
           obj.x_size = x_size_;
           obj.u_size = u_size_;
           obj.z_size = z_size_;
           
           % array init
           obj.F_array = zeros(x_size_, x_size_, h_size_);
           obj.H_array = zeros(z_size_, x_size_, h_size_);
           obj.y_tilde = zeros(z_size_, h_size_);
           obj.u_tilde = zeros(x_size_, h_size_);
           
           % function init
           obj.function_f = function_f_;
           obj.function_h = function_h_;
           obj.function_jf = function_jf_;
           obj.function_jh = function_jh_;
           
           % count init
           obj.count = 0;
           
           % init ok
           obj.is_init = "ok";
           r = obj.is_init;
       end
       
       function r = FIR_colson_run(obj,x_pre_, u_, y_)
           if obj.is_init == "ok"
               if obj.count == 0
                   argument_f = num2cell([0' u_']);
               else
                    argument_f = num2cell([x_pre_' u_']);
               end
               f_hat = obj.function_f(argument_f{:});
               argument_h = num2cell([f_hat']);
               h_hat = obj.function_h(argument_h{:});
               
               F =  obj.function_jf(argument_f{:});
               H =  obj.function_jh(argument_h{:});
               
               % accumulating array
               obj.F_array(:,:,1:obj.h_size - 1) = obj.F_array(:,:,2:obj.h_size);
               obj.F_array(:,:,obj.h_size) = F;
               obj.H_array(:,:,1:obj.h_size-1) = obj.H_array(:,:,2:obj.h_size);
               obj.H_array(:,:,obj.h_size) = H;
               obj.y_tilde(:,1:obj.h_size-1) = obj.y_tilde(:,2:obj.h_size);
               obj.y_tilde(:,obj.h_size) = y_ - (h_hat - H * f_hat);
               obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
               obj.u_tilde(:,obj.h_size) = f_hat - F * x_pre_;
               
               % calculate x_hat
               if obj.count > obj.h_size
               end
               obj.count = obj.count + 1;
           else
               disp("you must init class\r\nCall (filtering_init(obj, ...)\r\n");
           end
       end
       
   end
end