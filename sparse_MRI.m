function sparse_MRI
close all
%[X,map] = imread('circle_of_willis_tof_mip_s.gif');
[X, map] = imread('sparse_retinal_angiogram .gif');
[M,N,T] = size(X);
sparsity = nnz(X)/(M*N*T)
figure(1)
subplot(1,2,1)
imshow(X)
X_d = double(X);


% [~,threshold] = edge(X_d,'sobel');
% fudgeFactor = 0.5;
% X_d = edge(X_d,'sobel',threshold * fudgeFactor);
X_d = tensor(X_d);

T = tucker_als(X_d,[30 35]);
rec = uint8(double(T));
figure(1)
subplot(1,2,2)
imshow(rec)
%imshow(X)
end