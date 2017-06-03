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

% eigen decomposition of laplacian
[V, E] = eig(L) ;

% Further process is done only on second eigen vector for now
% It can be done for some k eigen vectors having small eigen values

% two halves of second eigen vector
% V(542, 2) = -0.0060 ;
E1 = V(1:nf1, 2) ;
E2 = V(nf1+1:nf, 2) ;

% first image interpolation
c1 = unique(f1(1, :)) ;
r1 = unique(f1(2, :)) ;

[C1, R1] = meshgrid(c1, r1) ;
ER1 = reshape(E1, [size(r1, 2) size(c1, 2)]) ;

[NC1, NR1] = meshgrid(1:n1, 1:m1) ;

I1 = interp2(C1, R1, ER1, NC1, NR1, 'linear', 0) ;

% second image interpolation
c2 = unique(f2(1, :)) ;
r2 = unique(f2(2, :)) ;

[C2, R2] = meshgrid(c2, r2) ;
ER2 = reshape(E2, [size(r2, 2) size(c2, 2)]) ;

[NC2, NR2] = meshgrid(1:n2, 1:m2) ;

I2 = interp2(C2, R2, ER2, NC2, NR2, 'linear', 0) ;


figure ;
subplot(1, 2, 1) ;
imagesc(I1) ;
subplot(1, 2, 2) ;
imagesc(I2) ;