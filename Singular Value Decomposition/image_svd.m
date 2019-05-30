close all
clear all
clc

%reading and converting the image
inImage=imread('image/iron_man.jpg');
inImage=rgb2gray(inImage);
inImageD=double(inImage);

myinImage=imread('image/iron_man.jpg');
myinImageD = double(myinImage);

%sort RGB
[nx, ny, n] = size(myinImageD)
result = zeros(nx,ny,n);
if n == 3
    R = myinImageD(:,:,1);
    G = myinImageD(:,:,2);
    B = myinImageD(:,:,3);
end

% decomposing the image using singular value decomposition
[U,S,V]=svd(inImageD);

[UR, SR, VR] = svd(R);
[UG, SG, VG] = svd(G);
[UB, SB, VB] = svd(B);


% Using different number of singular values (diagonal of S) to compress and
% reconstruct the image
dispEr = [];
numSVals = [];

dispEr2 = [];
numSVals2 = [];
    figure;
for N=50:30:60
    % store the singular values in a temporary var
    C = S;
    CR = SR;
    CB = SB;
    CG = SG;
    
    % discard the diagonal values not required for compression
    C(N+1:end,:)=0;
    C(:,N+1:end)=0;
    
    CR(N+1:end,:)=0;
    CR(:,N+1:end)=0;
    CG(N+1:end,:)=0;
    CG(:,N+1:end)=0;
    CB(N+1:end,:)=0;
    CB(:,N+1:end)=0;
    % Construct an Image using the selected singular values
    D=U*C*V';
    
    DR = UR * CR * VR';
    DG = UG * CG * VG';
    DB = UB * CB * VB';
    
    result(:,:,1) = DR;
    result(:,:,2) = DG;
    result(:,:,3) = DB;

    % display and compute error
%     figure;
%     buffer = sprintf('Image output using %d singular values', N)
%     imshow(uint8(D));
%     title(buffer);
    error=sum(sum((inImageD-D).^2));
    

    % store vals for display
    dispEr = [dispEr; error];
    numSVals = [numSVals; N];
    
    
    buffer2 = sprintf('Image output using %d singular values', N);
    imshow(uint8(result));
    title(buffer2);
    pause(0.5);
end

% dislay the error graph

% figure; 
% title('Error in compression');
% plot(numSVals, dispEr);
% grid on
% xlabel('Number of Singular Values used');
% ylabel('Error between compress and original image');