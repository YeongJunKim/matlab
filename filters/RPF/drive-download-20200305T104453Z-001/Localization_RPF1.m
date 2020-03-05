function [state_hat] = Localization_RPF(f,h,z,control_commands)
% Localization based on tag-anchor distances using RPF (수정중)
%%% Input arguments %%%
% -> f,F : State equations, Jacobian (Symbolic, 현재 변수 대입 필요)
% -> h,H : Measurement equations, Jacobian (Symbolic, 현재 변수 대입 필요)
% -> z : Tag 1로부터의 거리값들, [d12 d13 d14]'
% -> control_commands : [delta_d delta_theta]'

%%% Output arguments %%%
% -> state_hat(1) = x_coor : 위치 x 좌표값
% -> state_hat(2) = y_coor : 위치 y 좌표값
% -> state_hat(3) = heading : heading angle

% syms x y theta u1 u2

%% Parameters
persistent N_particle G G_order Q R firstRun
persistent dim_state dim_measurement
persistent x_hat_particle q
if isempty(firstRun)
	N_particle = 100;
    dim_state = 3;
    dim_measurement = 3;
    G = eye(dim_state); G_order = size(G);
    
    sigma_x = 0.1; sigma_y = 0.1; sigma_theta = pi/10;
    sigma_z1 = 0.2; sigma_z2 = 0.2; sigma_z3 = 0.2;
    Q = blkdiag(sigma_x^2,sigma_y^2,sigma_theta^2);
    R = blkdiag(sigma_z1^2,sigma_z2^2,sigma_z3^2);
	
    x_hat_particle = ones(dim_state, N_particle);   % 각 파티클의 추정값
    q = N_particle^(-1)*ones(1,N_particle);
    
    firstRun = 1;
end

%% RPF
for i = 1 : N_particle
    if(G_order(1,2) == 1)           % w is scalar
        noise_w = sqrt(Q)*G*randn;
        noise_v = sqrt(R)*randn(dim_measurement,1);
    else
        noise_w = sqrt(Q)*G*randn(dim_state,1);
        noise_v = sqrt(R)*randn(dim_measurement,1);
    end
	
    x_hat_particle(:,i) = double(f(x_hat_particle(1,i),x_hat_particle(2,i),x_hat_particle(3,i),control_commands(1),control_commands(2)))' + noise_w;
    z_hat = double(h(x_hat_particle(1,i),x_hat_particle(2,i),x_hat_particle(3,i),control_commands(1),control_commands(2)))' + noise_v;
        
    error = z - z_hat;
    q(i)=1 / ( (2*pi)^(dim_measurement/2)*det(R)^(1/2) )*exp( -error'/inv(R)*error/2 );
end

q_sum = sum(q);
q_normalize = q / q_sum; % normalize weight

%% Resampling
i=1;
q_cumsum=cumsum(q_normalize,2);
Num_random = N_particle^(-1)*rand; 
U = zeros(1,N_particle);
temp = x_hat_particle;

for j = 1 : N_particle
    U(1,j) = Num_random + N_particle^(-1)*(j-1);
    while U(1,j) > q_cumsum(1,i)
         i = i + 1;
    end
    x_hat_particle(:,j) = temp(:,i);
end

%% Regularization
h_opt = ((4/5)^(1/9)) * (N_particle^(-1/9));

S = cov(x_hat_particle(:,:)');
det_S = det(S);

if det_S < 10^(-20)
	D = zeros(dim_state,dim_state);
else
	D = sqrtm(S);
end

for j = 1 : N_particle
	e_x = randn;
	e_y = randn;
    e_yaw = randn;

	jitter = h_opt*D*[e_x e_y e_yaw]';
	x_hat_particle(:,j) = x_hat_particle(:,j) + jitter;
end

state_hat = mean( x_hat_particle, 2 );

end