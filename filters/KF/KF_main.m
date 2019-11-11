function [x_hat] = KF_main(A_, H_, Q_, R, x_init, P, x)
    persistent firstRun
    if isempty(firstRun)
       firstRun = 1;
       
    end
    
    
    xp = A*x;
    Pp = A*P*A'+Q;
    K = Pp*H'*inv(H*Pp*H'+R);
    x = xp+K*(z-H*xp);
    P = Pp-K*H*Pp;    
end