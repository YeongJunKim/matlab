clear all
A = [2/3 1/3;
    1/3 2/3];
 
    

A_power=A^1000

[V D W] = eig(A)


%% next e.g.,
A = [7/12 1/12 1/12 1/12 1/12 1/12;
    7/12 1/12 1/12 1/12 1/12 1/12;
    7/12 1/12 1/12 1/12 1/12 1/12;
    1/12 1/12 1/12 1/12 1/12 7/12;
    1/12 1/12 1/12 1/12 1/12 7/12;
    1/12 1/12 1/12 1/12 1/12 7/12];

x_zero = [1 0 0 0 0 -1]';

A*x_zero

A_power = A ^1000

[V D W] = eig(A)



%% circulant balancing
n = 15
k = 0.4

A = zeros(n,n);

for i=1:n
   for j=1:n
        if i==j
            A(i,j) = 1-2*k;
        end
        if (i+1) == j
            A(i,j) = k;
        end
        if i == (j+1)
            A(i,j) = k;
        end
   end
end
A(1,n) = k;
A(n,1) = k;

A_power = A;
A_power = A_power*A;
