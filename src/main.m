clear all ;

path1 = '1.jpg' ;
path2 = '2.jpg' ;

% extract dense sift features from both images
fprintf('dense sift features for %s\n', path1) ;
[f1, d1, m1, n1] = get_sift_features(path1) ;
fprintf('dense sift features for %s\n', path2) ;
[f2, d2, m2, n2] = get_sift_features(path2) ;

% number of feature points of image1 and image2 is nf1 and nf2 respectively
nf1 = size(d1, 2) ;
nf2 = size(d2, 2) ;
nf = nf1 + nf2 ;

% create adjacency matrix withing the vertices of each image and across
% the vertices of two images and making complete adjacency matrix
d1 = double(d1') ;
d2 = double(d2') ;
Adj1 = squareform(pdist(d1, 'cosine')) ;
Adj2 = squareform(pdist(d2, 'cosine')) ;
Adjcross = pdist2(d1, d2, 'cosine') ;

A_upper = horzcat(Adj1, Adjcross) ;
A_lower = horzcat(Adjcross', Adj2) ;
A = vertcat(A_upper, A_lower) ;

% create laplacian matrix using formule L = D - A
D = sum(A) ;
L = -A ;
L(1:nf+1:nf*nf) = D(:) ;

[V, E] = eig(L) ;

k1 = V(1:nf1, 2) ;
k2 = V(nf1+1:nf, 2) ;

[M1, M2] = meshgrid(1:n1, 1:m1) ;

