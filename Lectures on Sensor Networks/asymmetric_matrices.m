A = [1 0 0 0 0 0;
                1 0 0 0 0 0;
                1 0 0 0 0 0;
                0 0 0 0 0 1;
                0 0 0 0 0 1;
                0 0 0 0 0 1];
            
x_zero = [1 0 0 0 0 -1]';
A_power = A^1000

[V D W] = eig(A)


x_zero 
x_1 = A*x_zero


x_deposit = x_zero
x_deposit = [x_deposit, x_1]