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

A = adjacency_mat(d1, d2) ;
[V, E] = eigen_spectram(A) ;

% Further process is done only on second eigen vector for now
% It can be done for some k eigen vectors having small eigen values

% two halves of second eigen vector
num = 2 ;
E1 = V(1:nf1, num) ;
E2 = V(nf1+1:nf, num) ;

I1 = reconstruction(E1, f1, m1, n1) ;
I2 = reconstruction(E2, f2, m2, n2) ;

figure ;
subplot(1, 2, 1) ;
imagesc(I1) ;
subplot(1, 2, 2) ;
imagesc(I2) ;
