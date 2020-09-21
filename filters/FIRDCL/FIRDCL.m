% paper: 
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-09-08
% Finite Impulse Response Filtering Algorithm with Distributed Cooperative
% Localization
% i-th input      : u1 u2 u3  ...    u(nh)    u(nh+1)  ... u(nh+nh2)...
% i-th measurement:    y1 y2 y3  ... y(nh)    y(nh+1)  ... y(nh+nh2)...
% j-th x_hat                         x_hat(1) x_hat(2) ... x_hat(nh2)...
% i-th x_hathat                      x_hat(1) x_hat(2) ... x_hat(nh2)...             

classdef FIRDCL < handle
   properties 
       
       % functions
       function_f;
       function_jf;
       function_h;
       function_jh; 
       function_r;
       
       % measurement
       z = [];
       % control input
       u = [];
       % esitimation x
       x_hat = [];
       x_pre = [];
       
       % horizon size
       h_size;
       h2_size;
       % x size
       x_size;
       % u size
       u_size;
       % z size
       z_size;
       % r size
       r_size;
        
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
       
       % counting
       count = 0;
       
       is_init = "no";
       
       % data saving
       
       x_appended;
       
   end
   methods
       %% function area
       
       function obj = FIRDCL(nh1_, nh2_, nx_, nz_, nr_, nu_, ... % 5
                                 F_, J_F_, H_, J_H_, R_, ... % 5
                                init_state_, neighbor_num_)  % 2
           % size mapping
           obj.h_size = nh1_; % horizon size
           obj.h2_size = nh2_; % relative measurement horizon size
           obj.x_size = nx_; % state size
           obj.u_size = nu_; % input size
           obj.z_size = nz_; % measurement size
           obj.r_size = nr_; % relative measurement size
           
           % array init
           obj.F_array = zeros(nx_, nx_, nh1_);
           obj.H_array = zeros(nz_, nx_, nh1_);
           obj.y_tilde = zeros(nz_, nh1_);
           obj.u_tilde = zeros(nx_, nh1_);
           
           % function init
           obj.function_f = F_;
           obj.function_h = H_;
           obj.function_jf = J_F_;
           obj.function_jh = J_H_;
           obj.function_r = R_;
           
           % count init
           obj.count = 2;
           
           % state init
           obj.init_state = init_state_;
           obj.x_pre = init_state_;
           
           % data saving
           obj.x_appended = zeros(nx_, []);
           obj.x_appended(:,1) = obj.x_pre;
           
           % init ok
           obj.is_init = "ok";
       end
       
       function r = estimate(obj, u_, z_, r_)
           if obj.is_init ~= "ok"
                error("No initialized");
           end
           argument_f = num2cell([x_pre' u_']);
           f_hat = obj.function_f(argument_f{:});
           argument_h = num2cell([f_hat']);
           h_hat = obj.function_h(argument_h{:});
           
           F = obj.function_jf(argument_f{:});
           H = obj.function_jh(argument_h{:});
           
           % accumulating array
           obj.F_array(:,:,1:obj.h_size - 1) = obj.F_array(:,:,2:obj.h_size);
           obj.F_array(:,:,obj.h_size) = F;
           obj.H_array(:,:,1:obj.h_size-1) = obj.H_array(:,:,2:obj.h_size);
           obj.H_array(:,:,obj.h_size) = H;
           obj.y_tilde(:,1:obj.h_size-1) = obj.y_tilde(:,2:obj.h_size);
           obj.y_tilde(:,obj.h_size) = z_ - (h_hat - H * f_hat);
           obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
           obj.u_tilde(:,obj.h_size) = f_hat - F * x_pre_;
           obj.r_tilde(:,1:obj.h2_size-1) = obj.r_tilde(:,2:obj.h_size);
           
           
       end
       
       function r = FIR_run(obj,x_pre_, u_, z_)
           if obj.is_init == "ok"
               argument_f = num2cell([x_pre_' u_']);
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
               obj.y_tilde(:,obj.h_size) = z_ - (h_hat - H * f_hat);
               obj.u_tilde(:,1:obj.h_size-1) = obj.u_tilde(:,2:obj.h_size);
               obj.u_tilde(:,obj.h_size) = f_hat - F * x_pre_;
               
               % calculate x_hat
               if obj.count > obj.h_size
                    r = FIR_main(obj.F_array,obj.H_array,obj.y_tilde,obj.u_tilde,obj.h_size);
               else
                   r = obj.f_hat;
               end
               obj.x_appended(:,obj.count) = r;
               obj.count = obj.count + 1;
           else
               error("No initialized");
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