function [f, d, M, N] = get_sift_features2(im_path, max_pixels, step, bin_size)
% Resize image so that total number of feature points doesn't exceed much.
% There will be almost `max_pixels/(step*step)` feature points for small bin_size.
% 
% Extract SIFT features for uniformly distributed keypoints for two
% scales( two different bin_size) and concate them to make a single 256D
% feature.

im = imread(im_path) ;
if size(im, 3) == 3
    im = rgb2gray(im) ;
end

[M, N] = size(im) ;
im = imresize(double(im), sqrt(max_pixels / (M*N)) ) ;
im = im2single(im) ;
[M, N] = size(im) ;

% frames to compute sift features
[p, q] = meshgrid(1:step:N, 1:step:M) ;
f = [p(:) q(:)]' ;
f(3, :) = 3 ;
f(4, :) = 0 ;

% sift at two bin_sizes
[f, d1] = vl_sift(im, 'frames', f, 'WindowSize', bin_size(1), 'verbose') ;
[f, d2] = vl_sift(im, 'frames', f, 'WindowSize', bin_size(2), 'verbose') ;

f = f(1:2, :) ;
d = vertcat(d1, d2) ;

fprintf('total number of features : %d\n', size(d, 2)) ;

end