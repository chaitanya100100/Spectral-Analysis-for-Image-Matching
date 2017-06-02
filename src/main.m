clear all ;

path1 = '1.jpg' ;
path2 = '2.jpg' ;

% extract dense sift features from both images
fprintf('dense sift features for %s\n', path1) ;
[f1, d1] = get_sift_features(path1) ;
fprintf('dense sift features for %s\n', path2) ;
[f2, d2] = get_sift_features(path2) ;

% nf1 will be number of vertices in graph from image1
% nf2 will be number of vertices in graph from image2
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
A = -A ;
A(1:nf+1:nf*nf) = D(:) ;