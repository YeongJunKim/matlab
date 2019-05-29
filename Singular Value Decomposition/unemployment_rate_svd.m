yr = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]';
pop = [4.33, 4.20, 3.93, 4.41, 4.63, 4.78, 4.98 4.89, 4.51, 5.66, 5.62, 4.49, 4.25, 5.02, 5.52, 5.25, 5.96, 6.26, 6.42]';

A = [yr.^4 yr.^3 yr.^2 yr ones(length(yr),1)];

B = pop;

[U D V] = svd(A);

figure(1);
% pseudo inverse with full SVD
subplot(1,3,1);
D_inv = D;
D_inv(1,1) = 1/D(1,1)
D_inv(2,2) = 1/D(2,2)
D_inv(3,3) = 1/D(3,3)
A_pinv = V*D_inv'*U';
X = A_pinv*B;
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
pop_approx = A*X;
 plot(yr, pop, '*', yr, pop_approx);
