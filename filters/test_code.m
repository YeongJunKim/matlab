clear all

% addpath("C:\Users\colson\Desktop\workspace\2019_대한전기학회_KIEE");

addpath('EKF');
addpath('FIR');
addpath('TEMP');
addpath('sample data');
%% load data
load('190503_multirobot_data.mat');
num = size(data.time_exp,2);

%% Anchor positions, sampling time
a = 20;
x1 = 0; y1 = 0;
x2 = 0; y2 = a;
x3 = a; y3 = a;
x4 = a; y4 = 0;

dt = 0.5;   % Sampling time

%% State and Measurement equations
syms x y theta
syms v_l v_theta
syms x_sub y_sub theta_sub
syms v_l_sub v_theta_sub
% State of main robot
f_main.x(x,y,theta,v_l,v_theta) = x + v_l * cos(theta) * dt;
f_main.y(x,y,theta,v_l,v_theta) = y + v_l * sin(theta) * dt;
f_main.theta(x,y,theta,v_l,v_theta) = theta + v_theta * dt;
f_main = [f_main.x f_main.y f_main.theta]';
Jacobian_F_main = jacobian(f_main,[x,y,theta]);
f_main = matlabFunction(f_main);
Jacobian_F_main = matlabFunction(Jacobian_F_main);

% State of sub robots
f_sub.x_main(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = x;
f_sub.y_main(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = y;
f_sub.x_sub(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = x_sub + v_l_sub * cos(theta_sub) * dt;
f_sub.y_sub(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = y_sub + v_l_sub * sin(theta_sub) * dt;
f_sub.theta(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = theta_sub + v_theta_sub * dt;
f_sub = [f_sub.x_main f_sub.y_main f_sub.x_sub f_sub.y_sub f_sub.theta]';
Jacobian_F_sub = jacobian(f_sub,[x,y,x_sub,y_sub,theta_sub]);
f_sub = matlabFunction(f_sub);
Jacobian_F_sub = matlabFunction(Jacobian_F_sub);

% Measurement of main robot
h_main.dist1(x,y,theta,v_l,v_theta) = sqrt((x-x1)^2 + (y-y1)^2);
h_main.dist2(x,y,theta,v_l,v_theta) = sqrt((x-x2)^2 + (y-y2)^2);
h_main.dist3(x,y,theta,v_l,v_theta) = sqrt((x-x3)^2 + (y-y3)^2);
h_main.dist4(x,y,theta,v_l,v_theta) = sqrt((x-x4)^2 + (y-y4)^2);
h_main.theta(x,y,theta,v_l,v_theta) = theta;
h_main = [h_main.dist1 h_main.dist2 h_main.dist3 h_main.dist4 h_main.theta]';
Jacobian_H_main = jacobian(h_main,[x,y,theta]);
h_main = matlabFunction(h_main);
Jacobian_H_main = matlabFunction(Jacobian_H_main);

% Measurement of sub robots
h_sub.x_main(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = x;
h_sub.y_main(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = y;
h_sub.d_x(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = x - x_sub;
h_sub.d_y(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = y - y_sub;
h_sub.theta(x,y,x_sub,y_sub,theta_sub,v_l_sub,v_theta_sub) = theta_sub;
h_sub = [h_sub.x_main h_sub.y_main h_sub.d_x h_sub.d_y h_sub.theta]';
Jacobian_H_sub = jacobian(h_sub,[x,y,x_sub,y_sub,theta_sub]);
h_sub = matlabFunction(h_sub);
Jacobian_H_sub = matlabFunction(Jacobian_H_sub);








%% Localization
i = 1;

x_main_EKF_data = zeros(3,num);
x_main_PEFFME_data = zeros(3,num);
x_sub1_EKF_data = zeros(5,num);
x_sub1_PEFFME_data = zeros(5,num);
x_sub2_EKF_data = zeros(5,num);
x_sub2_PEFFME_data = zeros(5,num);

x_main_hat_EKF = [10 10 0]';
x_main_hat_PEFFME = [10 10 0]';
x_sub1_hat_EKF = [10 10 10 10 0]';
x_sub1_hat_PEFFME = [10 10 10 10 0]';
x_sub2_hat_EKF = [10 10 10 10 0]';
x_sub2_hat_PEFFME = [10 10 10 10 0]';

d_x_main_EKF_data = zeros(3,num);
d_x_main_PEFFME_data = zeros(3,num);
d_x_sub1_EKF_data = zeros(5,num);
d_x_sub1_PEFFME_data = zeros(5,num);
d_x_sub2_EKF_data = zeros(5,num);
d_x_sub2_PEFFME_data = zeros(5,num);

d_x_main_hat_EKF = [10 10 0]';
d_x_main_hat_PEFFME = [10 10 0]';
d_x_sub1_hat_EKF = [10 10 10 10 0]';
d_x_sub1_hat_PEFFME = [10 10 10 10 0]';
d_x_sub2_hat_EKF = [10 10 10 10 0]';
d_x_sub2_hat_PEFFME = [10 10 10 10 0]';

%% FIR filter class
FIR_FILTER = [];
for i = 1:3
    FIR_FILTER(i).filter = FIR;
end
for i = 1:3
    EKF_FILTER(i).filter = EKF; 
end
% argument
% obj,% class object
% h_size_,          % hrozion size
% x_size_,          % state_size
% z_size_,          % measurement size
% u_size_,          % control input size
% function_f_,      % function of F matrix
% function_jf_,     % jacobian of F
% function_h_,      % function of H matrix
% function_jh_,     % jacobian of H
FIR_init(FIR_FILTER(1).filter, 6, 3, 5, 3, f_main, Jacobian_F_main, h_main, Jacobian_H_main, x_main_hat_PEFFME);
FIR_init(FIR_FILTER(2).filter, 6, 5, 5, 3, f_sub, Jacobian_F_sub, h_sub, Jacobian_H_sub, x_sub1_hat_PEFFME);
FIR_init(FIR_FILTER(3).filter, 6, 5, 5, 3, f_sub, Jacobian_F_sub, h_sub, Jacobian_H_sub, x_sub2_hat_PEFFME);

EKF_init(EKF_FILTER(1).filter, eye(3), (blkdiag(0.1, 0.1, 0.1)), (blkdiag(0.2,0.2,0.2,0.2,0.1)), f_main, Jacobian_F_main,h_main, Jacobian_H_main);
EKF_init(EKF_FILTER(2).filter, eye(5), (0.1 * eye(5)), (0.1*eye(5)), f_sub, Jacobian_F_sub, h_sub, Jacobian_H_sub);
EKF_init(EKF_FILTER(3).filter, eye(5), (0.1 * eye(5)), (0.1*eye(5)), f_sub, Jacobian_F_sub, h_sub, Jacobian_H_sub);

t0 = 0;
t1 = t0 + 30;
t2 = t1 + 5;
t3 = t2 + 30;

time1 = zeros(num,1);
time2 = zeros(num,1);
while(1)
    disp(i);
    
    d1 = data.measurement_main(1,i);
    d2 = data.measurement_main(2,i);
    d3 = data.measurement_main(3,i);
    d4 = data.measurement_main(4,i);
    theta_main = data.measurement_main(5,i);
    
    if d1*d2*d3*d4 == 0
        alpha = 0;
    else
        alpha = 1;
    end
    
    z_main = [d1 d2 d3 d4 theta_main]';
    z_sub1_EKF = [x_main_hat_EKF(1:2)' data.measurement_sub1(:,i)']';
    z_sub1_PEFFME = [x_main_hat_PEFFME(1:2)' data.measurement_sub1(:,i)']';
    z_sub2_EKF = [x_main_hat_EKF(1:2)' data.measurement_sub2(:,i)']';
    z_sub2_PEFFME = [x_main_hat_PEFFME(1:2)' data.measurement_sub2(:,i)']';
	
    
    if data.time_exp < t1
        control = [0.3 0]';
    elseif data.time_exp < t2
        control = [0 -pi/10]';
    else
        control = [0.3 0]';
    end
    
    
    %% new method
    tic();
    x_main_hat_EKF = EKF_run(EKF_FILTER(1).filter, x_main_hat_EKF, control, z_main);
    x_main_hat_PEFFME = FIR_PEFFME_run(FIR_FILTER(1).filter, x_main_hat_PEFFME, control, z_main, alpha);
    x_sub1_hat_EKF = EKF_run(EKF_FILTER(2).filter, x_sub1_hat_EKF, control, z_sub1_EKF);
    x_sub1_hat_PEFFME = FIR_PEFFME_run(FIR_FILTER(2).filter, x_sub1_hat_PEFFME, control, z_sub1_PEFFME, alpha);
    x_sub2_hat_EKF = EKF_run(EKF_FILTER(3).filter, x_sub2_hat_EKF, control, z_sub2_EKF);
    x_sub2_hat_PEFFME = FIR_PEFFME_run(FIR_FILTER(3).filter, x_sub2_hat_PEFFME, control, z_sub2_PEFFME, alpha);
    time1(i) = toc();
    
    %% existing method
    tic();
    d_x_main_hat_EKF = Localization_main_EKF(f_main,Jacobian_F_main,h_main,Jacobian_H_main,x_main_hat_EKF,z_main,control);
    d_x_main_hat_PEFFME = Localization_main_PEFFME(f_main,Jacobian_F_main,h_main,Jacobian_H_main,x_main_hat_PEFFME,z_main,control,alpha);
    d_x_sub1_hat_EKF = Localization_sub1_EKF(f_sub,Jacobian_F_sub,h_sub,Jacobian_H_sub,x_sub1_hat_EKF,z_sub1_EKF,control);
    d_x_sub1_hat_PEFFME = Localization_sub1_PEFFME(f_sub,Jacobian_F_sub,h_sub,Jacobian_H_sub,x_sub1_hat_PEFFME,z_sub1_PEFFME,control,alpha);
    d_x_sub2_hat_EKF = Localization_sub2_EKF(f_sub,Jacobian_F_sub,h_sub,Jacobian_H_sub,x_sub2_hat_EKF,z_sub2_EKF,control);
    d_x_sub2_hat_PEFFME = Localization_sub2_PEFFME(f_sub,Jacobian_F_sub,h_sub,Jacobian_H_sub,x_sub2_hat_PEFFME,z_sub2_PEFFME,control,alpha);
    time2(i) = toc();
    
    x_main_EKF_data(:,i) = x_main_hat_EKF;
    x_main_PEFFME_data(:,i) = x_main_hat_PEFFME;
    x_sub1_EKF_data(:,i) = x_sub1_hat_EKF;
    x_sub1_PEFFME_data(:,i) = x_sub1_hat_PEFFME;
    x_sub2_EKF_data(:,i) = x_sub2_hat_EKF;
    x_sub2_PEFFME_data(:,i) = x_sub2_hat_PEFFME;
    
    d_x_main_EKF_data(:,i) = x_main_hat_EKF;
    d_x_main_PEFFME_data(:,i) = x_main_hat_PEFFME;
    d_x_sub1_EKF_data(:,i) = x_sub1_hat_EKF;
    d_x_sub1_PEFFME_data(:,i) = x_sub1_hat_PEFFME;
    d_x_sub2_EKF_data(:,i) = x_sub2_hat_EKF;
    d_x_sub2_PEFFME_data(:,i) = x_sub2_hat_PEFFME;
    if i == num
        break;
    end
    i = i+1;
end

%% Estimation error analysis
error_main_EKF_data = zeros(3,num);
error_main_PEFFME_data = zeros(3,num);
error_sub1_EKF_data = zeros(3,num);
error_sub1_PEFFME_data = zeros(3,num);
error_sub2_EKF_data = zeros(3,num);
error_sub2_PEFFME_data = zeros(3,num);

interval_of_interest = 15:num;

for i = interval_of_interest
    error_main_EKF_data(:,i) = data.state_main(:,i) - x_main_EKF_data(:,i);
    error_main_PEFFME_data(:,i) = data.state_main(:,i) - x_main_PEFFME_data(:,i);
    error_sub1_EKF_data(:,i) = data.state_sub1(:,i) - x_sub1_EKF_data(3:5,i);
    error_sub1_PEFFME_data(:,i) = data.state_sub1(:,i) - x_sub1_PEFFME_data(3:5,i);
    error_sub2_EKF_data(:,i) = data.state_sub2(:,i) - x_sub2_EKF_data(3:5,i);
    error_sub2_PEFFME_data(:,i) = data.state_sub2(:,i) - x_sub2_PEFFME_data(3:5,i);
end

% Position error
position_error_main_EKF_data = sqrt(error_main_EKF_data(1,interval_of_interest).^2 + error_main_EKF_data(2,interval_of_interest).^2);
position_error_main_PEFFME_data = sqrt(error_main_PEFFME_data(1,interval_of_interest).^2 + error_main_PEFFME_data(2,interval_of_interest).^2);
position_error_sub1_EKF_data = sqrt(error_sub1_EKF_data(1,interval_of_interest).^2 + error_sub1_EKF_data(2,interval_of_interest).^2);
position_error_sub1_PEFFME_data = sqrt(error_sub1_PEFFME_data(1,interval_of_interest).^2 + error_sub1_PEFFME_data(2,interval_of_interest).^2);
position_error_sub2_EKF_data = sqrt(error_sub2_EKF_data(1,interval_of_interest).^2 + error_sub2_EKF_data(2,interval_of_interest).^2);
position_error_sub2_PEFFME_data = sqrt(error_sub2_PEFFME_data(1,interval_of_interest).^2 + error_sub2_PEFFME_data(2,interval_of_interest).^2);

% Position RMSE
position_RMSE_main_EKF = sqrt(mean(error_main_EKF_data(1,interval_of_interest).^2 + error_main_EKF_data(2,interval_of_interest).^2));
position_RMSE_main_PEFFME = sqrt(mean(error_main_PEFFME_data(1,interval_of_interest).^2 + error_main_PEFFME_data(2,interval_of_interest).^2));
position_RMSE_sub1_EKF = sqrt(mean(error_sub1_EKF_data(1,interval_of_interest).^2 + error_sub1_EKF_data(2,interval_of_interest).^2));
position_RMSE_sub1_PEFFME = sqrt(mean(error_sub1_PEFFME_data(1,interval_of_interest).^2 + error_sub1_PEFFME_data(2,interval_of_interest).^2));
position_RMSE_sub2_EKF = sqrt(mean(error_sub2_EKF_data(1,interval_of_interest).^2 + error_sub2_EKF_data(2,interval_of_interest).^2));
position_RMSE_sub2_PEFFME = sqrt(mean(error_sub2_PEFFME_data(1,interval_of_interest).^2 + error_sub2_PEFFME_data(2,interval_of_interest).^2));

% Heading angle RMSE
angle_RMSE_main_EKF = sqrt(mean(error_main_EKF_data(3,interval_of_interest).^2));
angle_RMSE_main_PEFFME = sqrt(mean(error_main_PEFFME_data(3,interval_of_interest).^2));
angle_RMSE_sub1_EKF = sqrt(mean(error_sub1_EKF_data(3,interval_of_interest).^2));
angle_RMSE_sub1_PEFFME = sqrt(mean(error_sub1_PEFFME_data(3,interval_of_interest).^2));
angle_RMSE_sub2_EKF = sqrt(mean(error_sub2_EKF_data(3,interval_of_interest).^2));
angle_RMSE_sub2_PEFFME = sqrt(mean(error_sub2_PEFFME_data(3,interval_of_interest).^2));

%% Plot
figure(1);
subplot(3,1,1);
plot(data.time_exp(interval_of_interest),position_error_main_EKF_data,'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(data.time_exp(interval_of_interest),position_error_main_PEFFME_data,'d-','color',[0.9 0.3 0.3]);
% ylim([0,2]);
xlabel('(a)', 'FontSize', 13);
ylabel('Estimation error (m)', 'FontSize', 13);
legend({'EKF','Proposed'},'Location','southeast','NumColumns',2, 'Location', 'northeast', 'FontSize', 13);
subplot(3,1,2);
plot(data.time_exp(interval_of_interest),position_error_sub1_EKF_data,'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(data.time_exp(interval_of_interest),position_error_sub1_PEFFME_data,'d-','color',[0.9 0.3 0.3]);
% ylim([0,1]);
xlabel('(b)', 'FontSize', 13);
ylabel('Estimation error (m)', 'FontSize', 13);
legend({'EKF','Proposed'},'Location','southeast','NumColumns',2, 'Location', 'northeast', 'FontSize', 13);

subplot(3,1,3);
plot(data.time_exp(interval_of_interest),position_error_sub2_EKF_data,'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(data.time_exp(interval_of_interest),position_error_sub2_PEFFME_data,'d-','color',[0.9 0.3 0.3]);
xlabel('(c)', 'FontSize', 13);
ylabel('Estimation error (m)', 'FontSize', 13);
legend({'EKF','Proposed'},'Location','southeast','NumColumns',2, 'Location', 'northeast', 'FontSize', 13);


figure(2);
subplot(2,1,1);
plot(x_main_EKF_data(1,interval_of_interest),x_main_EKF_data(2,interval_of_interest),'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(x_main_PEFFME_data(1,interval_of_interest),x_main_PEFFME_data(2,interval_of_interest),'d-','color',[0.9 0.3 0.3]);
plot(data.state_main(1,interval_of_interest),data.state_main(2,interval_of_interest),'xk'); 
axis equal
xlim([0 16]); ylim([6 20]);
xlabel('x (m)', 'FontSize', 13);
ylabel('y (m)', 'FontSize', 13);
legend({'EKF','Proposed','Real'},'Location','southeast', 'FontSize', 13);

subplot(2,1,2);
plot(x_sub1_EKF_data(3,interval_of_interest),x_sub1_EKF_data(4,interval_of_interest),'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(x_sub1_PEFFME_data(3,interval_of_interest),x_sub1_PEFFME_data(4,interval_of_interest),'d-','color',[0.9 0.3 0.3]);
plot(data.state_sub1(1,interval_of_interest),data.state_sub1(2,interval_of_interest),'xk'); 
axis equal
xlim([0 16]); ylim([6 20]);
xlabel('x (m)', 'FontSize', 13);
ylabel('y (m)', 'FontSize', 13);
legend({'EKF','Proposed','Real'},'Location','southeast', 'FontSize', 13);





figure(3);
plot(x_main_EKF_data(1,interval_of_interest),x_main_EKF_data(2,interval_of_interest),'*-','color',[0.3 0.3 0.9]); hold on; grid on;
plot(x_main_PEFFME_data(1,interval_of_interest),x_main_PEFFME_data(2,interval_of_interest),'d-','color',[0.9 0.3 0.9]); hold on; grid on;
plot(data.state_main(1,interval_of_interest),data.state_main(2,interval_of_interest),'b-x'); hold on; grid on;
axis equal

plot(x_sub1_EKF_data(3,interval_of_interest),x_sub1_EKF_data(4,interval_of_interest),'*-','color',[0.3 0.6 0.9]); hold on; grid on;
plot(x_sub1_PEFFME_data(3,interval_of_interest),x_sub1_PEFFME_data(4,interval_of_interest),'d-','color',[0.9 0.6 0.6]); hold on; grid on;
plot(data.state_sub1(1,interval_of_interest),data.state_sub1(2,interval_of_interest),'c-x'); hold on; grid on;
axis equal

plot(x_sub2_EKF_data(3,interval_of_interest),x_sub2_EKF_data(4,interval_of_interest),'*-','color',[0.6 0.5 0.1]); hold on; grid on;
plot(x_sub2_PEFFME_data(3,interval_of_interest),x_sub2_PEFFME_data(4,interval_of_interest),'d-','color',[0.6 0.5 0.6]);
plot(data.state_sub2(1,interval_of_interest),data.state_sub2(2,interval_of_interest),'r-x'); 
axis equal
xlim([0 16]); ylim([6 20]);
xlabel('x (m)', 'FontSize', 13);
ylabel('y (m)', 'FontSize', 13);
legend({'A_{EKF}','A_{Proposed}','A_{Real}','B_{EKF}','B_{Proposed}','B_{Real}','C_{EKF}','C_{Proposed}','C_{Real}'},'Location','southwest', 'FontSize', 13);
