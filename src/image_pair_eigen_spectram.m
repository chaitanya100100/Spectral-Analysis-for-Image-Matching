function [I1, I2] = image_pair_eigen_spectram(path1, path2, max_pixels, step, bin_size, num)
% This function computes joint eigen spectram of two images.
% This function's parameters are two image paths, max size of image, sift parameters
% and number of eigen functions to be returned.

% extract dense sift features from both images
fprintf('dense sift features for %s\n', path1) ;
[f1, d1, m1, n1] = get_sift_features(path1, max_pixels, step, bin_size) ;
fprintf('dense sift features for %s\n', path2) ;
[f2, d2, m2, n2] = get_sift_features(path2, max_pixels, step, bin_size) ;

% number of feature points of image1 and image2 is nf1 and nf2 respectively
nf1 = size(d1, 2) ;
nf2 = size(d2, 2) ;
nf = nf1 + nf2 ;

% joint graph spectram
A = adjacency_mat(d1, d2) ;
[V, E] = graph_eigen_spectram(A, num) ;

I1 = {} ;
I2 = {} ;

% reconstruction of spectral images
for i = 2 : num+1
    E1 = V(1:nf1, i) ;
    E2 = V(nf1+1:nf, i) ;

    I1{i-1} = reconstruction(E1, f1, m1, n1) ;
    I2{i-1} = reconstruction(E2, f2, m2, n2) ;
end

end