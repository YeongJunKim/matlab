%% SVD population estimation

yr = [1930, 1940, 1949, 1960, 1970, 1980, 1990, 2000, 2010]';
pop = [2044, 2355, 2017, 2499, 3144, 3741, 4339, 4599, 4799]';

A = [yr.^2 yr ones(length(yr),1)];

B = pop;

[U D V] = svd(A);


figure(1);
% pseudo inverse with full SVD
subplot(1,3,1);
D_inv = D;
D_inv(1,1) = 1/D(1,1);
D_inv(2,2) = 1/D(2,2);
D_inv(3,3) = 1/D(3,3);
A_pinv = V*D_inv'*U';
X = A_pinv*B;
X
[2019^2 2019 1]*X
pop_approx = A*X;
plot(yr, pop, '*', yr, pop_approx);

% pseudo inverse with truncated SVD (t=2)
subplot(1,3,2);
D_inv = D;
D_inv(1,1) = 1/D(1,1);
D_inv(2,2) = 1/D(2,2);
D_inv(3,3) = 0;
A_pinv = V*D_inv'*U';
X = A_pinv*B;
X
pop_approx = A*X;
plot(yr, pop, '*', yr, pop_approx);

% pseudo inverse with truncated SVD (t=1)
subplot(1,3,3);
D_inv = D;
D_inv(1,1) = 1/D(1,1);
D_inv(2,2) = 0;
D_inv(3,3) = 0;
A_pinv = V*D_inv'*U';
X = A_pinv*B;
X
pop_approx = A*X;
 plot(yr, pop, '*', yr, pop_approx);
