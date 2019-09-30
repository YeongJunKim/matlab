function [x_hat] = EKF_main(f, h, jac_f, jac_h, x_pre, z, u, P, Q, R)

persistent P_ Q_ R_ firstRun P_I

if isempty(firstRun)
    P_ = P;
    Q_ = Q;
    R_ = R;
    
    firstRun = 1;
    
    P_I = eye(size(P_));
end
    
%% Substituting variables into Jacobian matrices
args = num2cell([x_pre' u']);
F = jac_f(args{:});
H = jac_h(args{:});


%% EKF algorithm
x_hat_temp = f(args{:});
P_ = F * P_ * F' + Q_;
K = P_ * H' / (H*P_*H' + R_);

Inno = z - h(args{:});
x_hat = x_hat_temp + K * Inno;
P_ = (P_I - K*H) * P_ * (P_I - K*H)' + K*R_*K';

end

