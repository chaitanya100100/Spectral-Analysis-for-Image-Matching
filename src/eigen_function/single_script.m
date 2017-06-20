% This script combines, in a single script, all tasks done by main.m,
% adjacency_mat.m, get_sift_features.m, graph_eigen_spectram.m,
% reconstruction.m.

clear
close all;

im1 = imread('../dataset/01.png');
im2 = imread('../dataset/02.png');

if size(im1,3) == 3
	im1 = rgb2gray(im1);
end
if size(im2,3) == 3
    im2 = rgb2gray(im2);
end
im1 = im2single(im1);
im2 = im2single(im2);
im1 = imresize(im1,0.3);
im2 = imresize(im2,0.3);

[p1,d1] = vl_dsift(im1,'step',4,'size',10);
[p2,d2] = vl_dsift(im2,'step',4,'size',10);
X = double([d1,d2]'+1);
N1 = size(d1,2);
N2 = size(d2,2);
N = N1 + N2;

W = pdist2(X,X,'cosine');
D = diag(sum(W));
D_cap = D^(-1/2);
L = eye(N) - D_cap * W * D_cap;

k = 5;
[U,lam] = eigs(L,k);
U1 = U(1:N1,:);
U2 = U(N1+1:N1+N2,:);

[X1,Y1] = meshgrid(1:size(im1,2),1:size(im1,1));
[X2,Y2] = meshgrid(1:size(im2,2),1:size(im2,1));
for i=2:k
	out1{i-1} = griddata(p1(1,:)',p1(2,:)',abs(U1(:,i)),X1,Y1);
	out2{i-1} = griddata(p2(1,:)',p2(2,:)',abs(U2(:,i)),X2,Y2);
end

for i=1:4
	figure
    subplot(1,2,1),imagesc(out1{i}),colormap jet
    title(strcat('eigenvector',int2str(i+1)));
    subplot(1,2,2),imagesc(out2{i}),colormap jet
end
