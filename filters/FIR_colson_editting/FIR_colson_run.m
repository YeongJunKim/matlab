
clear all


x = [];
x_hat = [];
y = [];
y_real = [];
cnt = 1;




%% 
f_main = [];
h_main = [];
cf_main = [];
vh_main = [];

x1 = 0; y1 = 0;
x2 = 0; y2 = 10;
x3 = 10;y3 = 10;
x4 = 10;y4 = 0;



%% main robot functions init

syms x y theta
syms v_l v_theta

dt = 0.01;


f_main.x(x,y,theta,v_l,v_theta) = x + v_l * cos(theta) * dt;
f_main.y(x,y,theta,v_l,v_theta) = y + v_l * cos(theta) * dt;
f_main_theta(x,y,theta,v_l,v_theta) = theta + v_theta * dt;
f_main = [f_main.x, f_main.y, f_main_theta]';
f_jacF = jacobian(f_main, [x,y,theta]);
f_main = matlabFunction(f_main);
f_jacF = matlabFunction(f_jacF);

h_main.dist1(x,y,theta,v_l,v_theta) = sqrt((x-x1)^2 + (y-y1)^2);
h_main.dist2(x,y,theta,v_l,v_theta) = sqrt((x-x2)^2 + (y-y2)^2);
h_main.dist3(x,y,theta,v_l,v_theta) = sqrt((x-x3)^2 + (y-y3)^2);
h_main.dist4(x,y,theta,v_l,v_theta) = sqrt((x-x4)^2 + (y-y4)^2);
h_main.theta(x,y,theta,v_l,v_theta) = theta;
h_main = [h_main.dist1, h_main.dist2, h_main.dist3, h_main.dist4, h_main.theta]';
h_jacF = jacobian(h_main,[x,y,theta]);
h_main = matlabFunction(h_main);
h_jacF = matlabFunction(h_jacF);


syms xx
syms yy
syms u
cf_main.y(xx, u) = cos(xx+u);
cf_main = [cf_main.y]';
cf_jacF = jacobian(cf_main, u);
cf_main = matlabFunction(cf_main);
cf_jacF = matlabFunction(cf_jacF);

vh_main.value(yy) = yy;
vh_main = [vh_main.value]';
% vh_jacF = jacobian(vh_main, yy);
vh_jacF = vh_main;
vh_main = matlabFunction(vh_main);
vh_jacF = matlabFunction(vh_jacF);


FIR_FILTER = [];
FIR_FILTER(1).filter = FIR_colson;
FIR_FILTER(2).filter = FIR_colson;
% FIR_FILTER(3).filter = EKF_colson;
% FIR_FILTER(40.filter = UKF_colson;
% FIR_colson_init(FIR_FILTER(1).filter, 6, 1, 1, 1, f_main,f_jacF,h_main,h_jacF);
FIR_colson_init(FIR_FILTER(2).filter, 6, 1, 1, 1, cf_main, cf_jacF, vh_main, vh_jacF);
% FIR_colson_init(FIR_FILTER(2).filter, 6, 1, 1, 1);

%% lcal robot functions init



%%
for i = 0:0.01:10
    disp(i);
    x(cnt) = i;
    [y(cnt) y_real(cnt)] = data_generating(x(cnt), "linear");
    
    
    
%     FIR_colson_run(FIR_FILTER(1).filter, 0.01, y(cnt));
    if i == 0
        FIR_colson_run(FIR_FILTER(2).filter,0, 0.01, y(cnt));
    else
        FIR_colson_run(FIR_FILTER(2).filter,x(cnt - 1), 0.01, y(cnt));
    end
%     x_hat = FIR_colson(x(cnt), y, 0.001);
    
     
    cnt = cnt + 1;
end
%%


figure(1);
subplot(1,2,1);
plot(x,y);
subplot(1,2,2);
plot(x,y_real);
