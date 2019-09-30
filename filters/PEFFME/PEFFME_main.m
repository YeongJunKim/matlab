function [x_hat] = PEFFME_main(f, h, jac_f, jac_h, x_pre, z, u, alpha, h_size)
    
persistent firstRun
persistent M F_array H_array z_array u_array y_tilde_array u_tilde_array
persistent dim_state dim_z
persistent state_tilde count

if isempty(firstRun)
    M = h_size;
    
    dim_state = size(x_pre, 1);
    dim_z = size(z,1);
    
    F_array = zeros(dim_state, dim_state, M);
    H_array = ones(dim_z, dim_state, M);
    z_array = zeros(dim_z, M);
    y_tilde_array = zeros(dim_z, M);
    u_array = zeros(length(u), M);
    u_tilde_array = zeros(dim_state, M);
    
    count = 0;
    state_tilde = x_pre;
    
    firstRun = 1;
end
%% Substituting variables into Jacobian matrices
argsf = num2cell([x_pre' u']);
F = jac_f(argsf{:});
H = jac_h(argsf{:});

f_hat = f(argsf{:});

argsh = num2cell([f_hat' u']);
h_hat = h(argsh{:});

%% Predict measurement if alpha = 0
if alpha == 0
    z = h_hat;
end

%% Array matrices
F_array(:,:,1:M-1) = F_array(:,:,2:M);
F_array(:,:,M) = F;
H_array(:,:,1:M-1) = H_array(:,:,2:M);
H_array(:,:,M) = H;
z_array(:,1:M-1) = z_array(:,2:M);
z_array(:,M) = z;
u_array(:,1:M-1) = u_array(:,2:M);
u_array(:,M) = u;
y_tilde_array(:,1:M-1) = y_tilde_array(:,2:M);
y_tilde_array(:,M) = z - (h_hat - H * f_hat);
u_tilde_array(:,1:M-1) = u_tilde_array(:,2:M);
u_tilde_array(:,M) = f_hat - F * x_pre;

count = count + 1;

%% result
if count > M
    x_hat = PEFFME(F_array, H_array, y_tilde_array, u_tilde_array, M);
else
    x_hat = state_tilde;
end

