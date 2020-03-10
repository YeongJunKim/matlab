clear all
clf
clc


%% add path & data
addpath('../sample data');
load('190503_multirobot_data.mat');
num = size(data.time_exp,2);

%% anchor positions, sampling time
a = 20;
x1 = 0; y1 = 0;
x2 = 0; y2 = a;
x3 = a; y3 = a;
x4 = a; y4 = 0;

dt = 0.5;   % Sampling time

%% Process & Observation equations
syms x y theta
syms v_l v_theta
syms x_sub y_sub theta_sub
syms v_l_sub v_theta_sub
syms u
syms v
syms k yk xk
% State of main robot
f_main.x(x,y,theta,v_l,v_theta) = x + v_l * cos(theta) * dt;
f_main.y(x,y,theta,v_l,v_theta) = y + v_l * sin(theta) * dt;
f_main.theta(x,y,theta,v_l,v_theta) = theta + v_theta * dt;
f_main = [f_main.x f_main.y f_main.theta]';
f_main = matlabFunction(f_main);
% Measurement of main robot
h_main.dist1(x,y,theta,v_l,v_theta) = sqrt((x-x1)^2 + (y-y1)^2);
h_main.dist2(x,y,theta,v_l,v_theta) = sqrt((x-x2)^2 + (y-y2)^2);
h_main.dist3(x,y,theta,v_l,v_theta) = sqrt((x-x3)^2 + (y-y3)^2);
h_main.dist4(x,y,theta,v_l,v_theta) = sqrt((x-x4)^2 + (y-y4)^2);
h_main.theta(x,y,theta,v_l,v_theta) = theta;
h_main = [h_main.dist1 h_main.dist2 h_main.dist3 h_main.dist4 h_main.theta]';
h_main = matlabFunction(h_main);
%% init data
x_main_PF_data = zeros(3,num);
x_main_PF_init = [3 18 0]';
x_main_PF_data(:,1) = x_main_PF_init;


samples = 200;
appended_num = 1000;
P = blkdiag(0.01, 0.01, 0.01);
Q = blkdiag(0.01, 0.01, 0.01);
R = blkdiag(0.2, 0.2, 0.2, 0.2, 0.1);
resampling_strategy = 'multinomial_resampling';


t0 = 0;
t1 = t0 + 30;
t2 = t1 + 5;
t3 = t2 + 30;

%% filter init
PF_FILTER = PF;

PF_init(PF_FILTER, samples, x_main_PF_data(:,1), appended_num, P, Q, R, f_main, h_main, resampling_strategy);

% start from i = 2 
% because init states are stored above "PF_init" function.
i = 2;
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
	
    if data.time_exp(i) < t1
        control = [0.3 0]';
    elseif data.time_exp(i) < t2
        control = [0 -pi/10]';
    else
        control = [0.3 0]';
    end
    
    x_main_PF_data(:,i) = PF_run(PF_FILTER, control, z_main);
    
    i = i + 1;
    if i == num
        break;
    end
end
    
interval_of_interest = 1:num-1;
plot(PF_FILTER.x_appended(1,interval_of_interest), PF_FILTER.x_appended(2,interval_of_interest)); hold on; grid on;
plot(data.state_main(1,interval_of_interest),data.state_main(2,interval_of_interest),'xk'); 
xlim([0,15])
ylim([5,20])


































