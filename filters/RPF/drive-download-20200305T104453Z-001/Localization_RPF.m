function [state_hat] = Localization_RPF(x_init,f,h,z,control_commands)
% Localization based on tag-anchor distances using RPF (수정중)
%%% Input arguments %%%
% -> f,F : State equations, Jacobian (Symbolic, 현재 변수 대입 필요)
% -> h,H : Measurement equations, Jacobian (Symbolic, 현재 변수 대입 필요)
% -> z : Tag 1로부터의 거리값들 및 heading angle, [d12 d13 d14 theta]'
% -> control_commands : [delta_d delta_theta]'

%%% Output arguments %%%
% -> state_hat(1) = x_coor : 위치 x 좌표값
% -> state_hat(2) = y_coor : 위치 y 좌표값
% -> state_hat(3) = heading : heading angle

%% Parameters
persistent N_particle G Q R firstRun
persistent sqrt_Q sqrt_R
persistent dim_state dim_z
persistent x_hat_particle q Particles v_n h_opt
if isempty(firstRun)
	N_particle = 50;
    dim_state = 3;
    dim_z = size(z,1);
    G = gpuArray(eye(dim_state));
    
%     sigma_x = 0.3; sigma_y = 0.3; sigma_theta = 0.3;
%     sigma_z1 = 0.5; sigma_z2 = 0.5; sigma_z3 = 0.5; sigma_z4 = 1;
%     Q = gpuArray(blkdiag(sigma_x^2,sigma_y^2,sigma_theta^2)); sqrt_Q = sqrt(Q);
    Q = gpuArray(blkdiag(0.01,0.01,0.01)); sqrt_Q = sqrt(Q);
%     R = gpuArray(blkdiag(sigma_z1^2,sigma_z2^2,sigma_z3^2,sigma_z4^2)); sqrt_R = sqrt(R);
    R = gpuArray(blkdiag(1.2,1.2,1.2,0.7)); sqrt_R = sqrt(R);
	
%     x_hat_particle = gpuArray(ones(dim_state, N_particle));   % 각 파티클의 추정값
    x_hat_particle = gpuArray(repmat(x_init,1,N_particle));
    q = gpuArray(N_particle^(-1)*ones(1,N_particle));
    for i = 1:N_particle
        Particles{i} = x_hat_particle(:,i);
    end
    
    v_n = 4*pi/3;
    h_opt = (8/v_n*(dim_state+4)*(2*sqrt(pi))^dim_state * N_particle)^(1/(dim_state+4));
    
    firstRun = 1;
end

%% RPF
u1 = control_commands(1);
u2 = control_commands(2);

Particles = cellfun(@(x) f(x(1),x(2),x(3),u1,u2) + sqrt_Q*G*randn(dim_state,1), Particles,'UniformOutput',false);
z_hat = cellfun(@(x) h(x(1),x(2),x(3),u1,u2) + sqrt_R*randn(dim_z,1), Particles,'UniformOutput',false);

for i=1:N_particle
    error = z - z_hat{i};
    q(i)=1 / ( (2*pi)^(dim_z/2)*det(R)^(1/2) )*exp( -error'/inv(R)*error/2 );        
end
q_sum = sum(q);
q = q / q_sum; % normalize weight
for i = 1:N_particle
	x_hat_particle(:,i) = Particles{i};
end

%% Resampling
N_eff = 1/sum(q.^2);
if N_eff < 0.1 * N_particle
    i=1;
    q_cumsum=cumsum(q,2);
    Num_random = N_particle^(-1)*rand; 
    U = zeros(1,N_particle);
    temp = Particles;

    for j = 1 : N_particle
        U(1,j) = Num_random + N_particle^(-1)*(j-1);
        while U(1,j) > q_cumsum(1,i)
             i = i + 1;
        end
        Particles{j} = temp{i};
    end

    for i = 1:N_particle
        x_hat_particle(:,i) = Particles{i};
    end

    q = gpuArray(N_particle^(-1)*ones(1,N_particle));
    
    %%% Regularization %%%
    S = cov(x_hat_particle');
    det_S = det(S);

    if det_S < 10^(-20)
        D = gpuArray.zeros(dim_state,dim_state);
    else
        D = chol(S,'lower');
    end

	for j = 1 : N_particle
        mu = mean(x_hat_particle,2);
        jitter = 1/(h_opt^dim_state) * Epanechnikov(mu/h_opt,v_n,dim_state) * randn(3,1);
        x_hat_particle(:,j) = x_hat_particle(:,j) + h_opt * D * jitter;
	end

    for i = 1:N_particle
        Particles{i} = x_hat_particle(:,i);
    end
end

%% Result

state_hat = gather(mean( x_hat_particle, 2 ));

end

function K = Epanechnikov(x,v_n,dim_state)
if norm(x) < 1
    K = (dim_state + 2)/(2*v_n)*(1 - norm(x)^2);
else
    K = 0;
end
end