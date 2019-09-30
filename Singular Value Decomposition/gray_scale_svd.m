clear all
%%
col = 10
row = 8
deletesize = 2
recsize = 10

if col == row
    rec = col
end
if col < row
    rec = col
end
if col > row
    rec = row
end


image = zeros(col * recsize, row * recsize);
image2 = zeros(col * recsize, row * recsize);

A = randi([10 100],col,row)


%%



[U D V] = svd(A)

D_inv = D

for i=0:deletesize
    for j=0:deletesize
        D_inv(rec-i,rec-j)=0;
    end
end

A_ = U*D_inv*V';

% original image

for i = 0:col-1
    for j = 0:row-1
        for k = 1:recsize
            for z = 1:recsize
                image(i*recsize+k,j*recsize+z) = A(i+1,j+1);
                image2(i*recsize+k,j*recsize+z) = A_(i+1,j+1);
            end
        end
    end
end


figure(1);
subplot(1,2,1);
imshow(uint8(image));

subplot(1,2,2);
imshow(uint8(image2));
% achieved image


