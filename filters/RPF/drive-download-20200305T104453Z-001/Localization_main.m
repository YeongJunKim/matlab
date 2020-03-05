%% Initializing
clear all
load('0706measurement_data.mat')
x_hat = [1 6 0.1]';

%% Position information of Anchors
distance = 0.6 * 10;
x1 = 0; y1 = 0;
x2 = 0; y2 = distance;
x3 = distance; y3 = distance;
x4 = distance; y4 = 0;

%% Construct Jacobian Matrices
syms x y theta      % state variables
syms u1 u2          % input variables, u1: 이동 거리, u2: 각도 변화

f1(x,y,theta,u1,u2) = x + u1 * cos(theta + 1/2 * u2);
f2(x,y,theta,u1,u2) = y + u1 * sin(theta + 1/2 * u2);
f3(x,y,theta,u1,u2) = theta + u2;
f = [f1 f2 f3]';
clear f1 f2 f3
Jacobian_F = jacobian(f,[x y theta]);

h1(x,y,theta,u1,u2) = sqrt((x-x1)^2 + (y-y1)^2) - sqrt((x-x2)^2 + (y-y2)^2);
h2(x,y,theta,u1,u2) = sqrt((x-x1)^2 + (y-y1)^2) - sqrt((x-x3)^2 + (y-y3)^2);
h3(x,y,theta,u1,u2) = sqrt((x-x1)^2 + (y-y1)^2) - sqrt((x-x4)^2 + (y-y4)^2);
h4(x,y,theta,u1,u2) = theta;
h = [h1 h2 h3 h4]';
clear h1 h2 h3 h4
Jacobian_H = jacobian(h,[x y theta]);

f = matlabFunction(f);
Jacobian_F = matlabFunction(Jacobian_F);
h = matlabFunction(h);
Jacobian_H = matlabFunction(Jacobian_H);

%% For Saving data
x_hat_data = zeros(4,100);
i = 1;
tic

%% Localization
while(1)
    %%% 센서 정보(d1,d2,d3,d4), 모바일 로봇으로의 control command 정보(delta_d, delta_theta)가 들어올 자리 %%%
%     z1 = d1 - d2;
%     z2 = d1 - d3;
%     z3 = d1 - d4;
%     z4 = theta;
%     z = [z1 z2 z3 z4]';
%     control = [delta_d delta_theta]';

    d1 = measurement_data(1,i);
    d2 = measurement_data(2,i);
    d3 = measurement_data(3,i);
    d4 = measurement_data(4,i);
    
    z1 = d1 - d2;
    z2 = d1 - d3;
    z3 = d1 - d4;
    z4 = 0;
    
    z = [z1 z2 z3 z4]';
    
%     z = [0 0 0 0]';
    control = [0 0]';
    
%     tic
	x_hat = Localization_RPF(x_hat,f,h,z,control)
%     toc
    
    x_hat_data(1,i) = i;
    x_hat_data(2:4,i) = x_hat(1:3);    
    
    if i == 100
        break;
    end
    i = i+1;
        
    %%% Plotting 또는 저장 %%%
    
end

%% Covariance 분석
cov(x_hat_data(2:4,:)')

%% Plot
figure(); plot(x_hat_data(1,:),x_hat_data(2,:)); ylim([0 distance]); title('x');
figure(); plot(x_hat_data(1,:),x_hat_data(3,:)); ylim([0 distance]); title('y');