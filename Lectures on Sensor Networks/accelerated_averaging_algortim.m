A = [1/2 1/2;
    1/2 1/2];

beta = 3
residual = 1-beta
A_ = eye(2) * beta * A


T_beta = [A_ (eye(2)*residual)
         eye(2) eye(2)*0];
     
T_beta_power = T_beta^1000
A_power = A^10000




A = [0 0 -2;
    0 1 0;
    1 0 3];

eig(A)