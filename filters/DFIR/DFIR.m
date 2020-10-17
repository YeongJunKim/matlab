% paper: 
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-07-28
% Distributed Finite Impulse Response Algorithm

classdef DFIR < handle
   properties 
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       function_r;
       function_jr; 
       
       % measurement
       z = [];
       % control input
       u = [];
       % esitimation x
       x_hat = [];
       x_pre = [];
       
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
       
       
       % init state
       init_state;
       
       % z_
       z_;
       % u_
       u_;
       
       % neighbors
       nn;
       
       % counting
       count = 0;
       
       is_init = "no";
       
       % data saving
       x_appended;
       % square error
       x_se;
       x_rmse;
   end
   methods
       %% function area
       
       % arguments
       % h_size_, x_size_, z_size, u_size_: size info.
       % function*: function of dynamics, measurement, relative measurement
       % nn: neighbor number
       function obj = DFIR(h_size_, x_size_, z_size_, u_size_, function_f_, function_jf_, function_h_, function_jh_, function_r_, function_jr_, init_state_, nn_)
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
           obj.count = 1;
           
           % state init
           obj.init_state = init_state_;
           obj.x_pre = init_state_;
           
           % data saving
           obj.x_appended = zeros(x_size_, []);
           
           % init ok
           obj.is_init = "ok";
       end
       
       
       % This function for the relative measurement 
       % u_ -> conrol input
       % pj_ -> row vector consist of neighbor's x and y, augmented.
       % z_ -> relative measurement, [x_(j=1), y_(j=1), x_(j=2),
       % y_(j=2),...x_(j=N_i), y_(j=N_i), theta]'
       
       function r = estimate2(obj, u_, pj_, z_)
           if strcmp(obj.is_init, "ok") == 0
               error("you must init class    : Call (filtering_init(obj, ...)");
           end
           argument_f = num2cell([obj.x_pre' u_']);
           f_hat = obj.function_f(argument_f{:});
           argument_h = num2cell([f_hat' pj_']);
           h_hat = obj.function_h(argument_h{:});
%            h_hat(4:7) = wrapTo2Pi(h_hat(4:7));
           
           F =  obj.function_jf(argument_f{:});
           H =  obj.function_jh(argument_h{:});
           
           % accumulating array
           obj.F_array(:,:,1:obj.h_size - 1) = obj.F_array(:,:,2:obj.h_size);
           obj.F_array(:,:,obj.h_size) = F;
           obj.H_array(:,:,1:obj.h_size-1) = obj.H_array(:,:,2:obj.h_size);
           obj.H_array(:,:,obj.h_size) = H;
           obj.y_tilde(:,1:obj.h_size-1) = obj.y_tilde(:,2:obj.h_size);
           obj.y_tilde(:,obj.h_size) = z_ - (h_hat - H * f_hat);
           obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
           obj.u_tilde(:,obj.h_size) = f_hat - F * obj.x_pre;
           
           % calculate x_hat
           if obj.count > obj.h_size
               r = FIR_main(obj.F_array,obj.H_array,obj.y_tilde,obj.u_tilde,obj.h_size);
           else
               r = f_hat;
           end
           obj.x_pre = r;
           obj.x_appended(:,obj.count) = r;
           obj.count = obj.count + 1;
       end
       
       function r = estimate(obj, u_, z_, pj_, alpha_)
           if obj.is_init == "ok"
%                disp("debug");
%                disp(obj.x_pre);

               argument_f = num2cell([obj.x_pre' u_']);
               f_hat = obj.function_f(argument_f{:});
               argument_h = num2cell([f_hat' pj_']);
               h_hat = obj.function_h(argument_h{:});
               % There is a matlab problem with [-Pi, Pi], so need to be
               % convert radian value go to [0, 2Pi].
               h_hat(3:4) = wrapTo2Pi(h_hat(3:4));
               
               F =  obj.function_jf(argument_f{:});
               H =  obj.function_jh(argument_h{:});
               z = (1 - alpha_) * h_hat + alpha_ * z_;
               
               % accumulating array
               obj.F_array(:,:,1:obj.h_size - 1) = obj.F_array(:,:,2:obj.h_size);
               obj.F_array(:,:,obj.h_size) = F;
               obj.H_array(:,:,1:obj.h_size-1) = obj.H_array(:,:,2:obj.h_size);
               obj.H_array(:,:,obj.h_size) = H;
               obj.y_tilde(:,1:obj.h_size-1) = obj.y_tilde(:,2:obj.h_size);
               obj.y_tilde(:,obj.h_size) = z - (h_hat - H * f_hat);
               obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
               obj.u_tilde(:,obj.h_size) = f_hat - F * obj.x_pre;
              
               % calculate x_hat
               if obj.count > obj.h_size
                    r = FIR_main(obj.F_array,obj.H_array,obj.y_tilde,obj.u_tilde,obj.h_size);
               else
                   r = f_hat;
               end
               obj.x_pre = r;
               obj.x_appended(:,obj.count) = r;
               obj.count = obj.count + 1;
           else
               error("you must init class    : Call (filtering_init(obj, ...)");
           end
       end
   end
end























function state_hat = FIR_main(F_array,H_array,y_tilde_array,u_tilde_array,M)
    A_big = Big_A(F_array,H_array,M);
    B_big = Big_B(F_array,H_array,M);
    C_big = Big_C(F_array,M);
    
    F0 = F_function(F_array,0,1,M);
    
    L = F0 / (A_big' * A_big) * A_big';
    M = -L * B_big + C_big;
    
   
    
    state_hat = L * reshape(y_tilde_array,[],1) + M * reshape(u_tilde_array,[],1);
end

function output = Big_A(F_array,H_array,M)
    m = size(H_array(:,:,1),1);
    n = size(H_array(:,:,1),2);
    
    output = zeros(M*m,n);
    
    for i = 1:M
        output(m*(i-1)+1:m*i,:) = Gamma_function(F_array,H_array,i,1,1);
    end
end

function output = Big_B(F_array,H_array,M)
    m = size(H_array(:,:,1),1);
    n = size(H_array(:,:,1),2);
    
    output = zeros(M*m,M*n);
    
    for i = 1:M-1
        for j = 1:M-1
            row_interval = m*i+1:m*(i+1);
            col_interval = n*(j-1)+1:n*j;
            
            output(row_interval,col_interval) = Gamma_function(F_array,H_array,i,j,0);
        end
    end
end

function output = Big_C(F_array,M)
    n = length(F_array(:,:,1));
    output = zeros(n,M*n);
    
    for i = 1:M
        output(:,n*(i-1)+1:n*i) = F_function(F_array,i,1,M);
    end
end

function output = Gamma_function(F_array,H_array,a,b,c)
    if a > b
        F_tmp = F_array(:,:,a-c);
        for i = 2:a-b
            F_tmp = F_tmp * F_array(:,:,a+1-c-i);
        end
        output = H_array(:,:,a+1-c) * F_tmp;
    elseif a == b
        output = H_array(:,:,a+1-c);
    else
        output = zeros(size(H_array(:,:,1)));
    end
end

function output = F_function(F_array,a,c,M)
    n = length(F_array(:,:,1));
    output = eye(n);
    
    if a < M
        for i = 1:M-a
            output = F_array(:,:,a+1+i-c) * output;
        end
    end
end