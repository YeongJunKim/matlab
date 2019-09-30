function [state_hat] = Localization_EKF(f,Jacobian_F,h,Jacobian_H,state_previous,z,control_commands)
% Localization based on tag-anchor distances using EKF
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
persistent P Q R firstRun
if isempty(firstRun)
	P = eye(3);
%     sigma_x = 0.3; sigma_y = 0.3; sigma_theta = 0.3;
%     sigma_z1 = 0.5; sigma_z2 = 0.5; sigma_z3 = 0.5; sigma_z4 = 1;
%     Q = blkdiag(sigma_x^2,sigma_y^2,sigma_theta^2);
    Q = blkdiag(0.01,0.01,0.01);
%     R = blkdiag(sigma_z1^2,sigma_z2^2,sigma_z3^2,sigma_z4^2);
    R = blkdiag(0.2,0.2,0.2,0.1);
    
    firstRun = 1;
end

%% Substituting variables into Jacobian matrices
F = Jacobian_F(state_previous(1),state_previous(2),state_previous(3),control_commands(1),control_commands(2));
H = Jacobian_H(state_previous(1),state_previous(2),state_previous(3),control_commands(1),control_commands(2));

%% EKF
% state_hat_temp = double(subs(f,[x,y,theta,u1,u2],[state_previous' control_commands']))';    % Prediction
state_hat_temp = f(state_previous(1),state_previous(2),state_previous(3),control_commands(1),control_commands(2));    % Prediction
P = F * P * F' + Q;
K = P * H' / (H*P*H' + R);
Inno = z - h(state_previous(1),state_previous(2),state_previous(3),control_commands(1),control_commands(2));  % Innovation
state_hat = state_hat_temp + K * Inno;   % Correction
P = (eye(3) - K*H) * P * (eye(3) - K*H)' + K*R*K';

end