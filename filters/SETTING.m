function fct = SETTING(select)

    addpath('PEFFME');
    addpath('EKF');
    addpath('KF');
    addpath('FIR');
    
    if select == "PEFFME"
        %% arguments
        % function [x_hat] = PEFFME_main(f, h, jac_f, jac_h, x_pre, z, u, alpha, h_size)
        % @f = function of f
        % @h = function of h
        % @jac_f = jacobian of function f
        % @jac_h = jacobian of function h
        % @x_pre = previous x_hat
        % @z = measurement data
        % @u = control input
        % @alpha = missing measurement flag (1 : missing) (0 : the other)
        % h_size = horizon size of FIR filter
        fct = @PEFFME_main;
    
    elseif select == "EKF"
        %% arguments
        % function [x_hat] = EKF_main(f, h, jac_f, jac_h, x_pre, z, u, P, Q, R) 
        % @f function of f
        % @h function of h
        % @jac_f = jacobian of function f
        % @jac_h = jacobian of function h
        % @x_pre = previous x_hat
        % @z = measurement data
        % @u = control input
        % @P = covariance
        % @Q,R = noise matrix
        
        fct = @EKF_main;
    elseif select == "KF"
        
        fct = @KF_main;
    end
        
end

