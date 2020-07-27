classdef DKFCL < handle
   properties 
       
       % algorithm type
       algorithm_num
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       
       A;
       B;
       C;
       D;
       
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
       filtering_type = "";
       first_run = 0;
   end
   methods
       %% function area
       function obj = DKFCL(P_, Q_, R_, A_, B_, C_, D_, init_, algorithm_num_)
           obj.algorithm_num = algorithm_num_;
           obj.P = P_;
           obj.Q = Q_;
           obj.R = R_;
           obj.A = A_;
           obj.B = B_;
           obj.C = C_;
           obj.D = D_;
           obj.x_appended = zeros(size(init_,1), []);
           obj.x_appended(:,1) = init_(:);
           obj.x_pre = init_;
           obj.count = 2;
           if(obj.algorithm_num ~= 1 && obj.algorithm_num ~= 2)
               error("Algorithm number can be the number as 1 or 2.");
           end
           obj.is_init = "ok";
       end
             
       function r = DKFCL_run(obj, u_, y_, z_)
           % paper title: Distributed Kalman filter for cooperative localization with integrated measurements
           % There are two algorithm: algorithm 1 focused on absoluted measurement is more important than relative measurement, algorithm 2 is versa.
           % u_ control input
           % y_ absolute measurement
           % z_ relative measurement
          if obj.is_init == "ok"
              % algorithm 1
              
              
              % algorithm 2
              
              
              obj.filtering_type = "DEKFCL";
          else
               error("you must init class    : Call (filtering_init(obj, ...)");
          end
       end
       
   end
end







