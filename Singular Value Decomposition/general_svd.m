clear all

A = [1 2;
    2 2];

[U D V] = svd(A)

U * D * V'

sqrt(5)