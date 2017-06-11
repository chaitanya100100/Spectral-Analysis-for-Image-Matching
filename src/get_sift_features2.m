function [f, d, M, N] = get_sift_features2(im_path, max_pixels, step, bin_size)
% Resize image so that total number of feature points doesn't exceed much.
% There will be almost `max_pixels/(step*step)` feature points for small bin_size.
% Extract image SIFT features using vlfeat.

%{
max_pixels = 48000 ;
step = 6 ;
bin_size = 6 ;
%}

im_raw = rgb2gray(imread(im_path)) ;
[M, N] = size(im_raw) ;
fprintf('image original size : [%d, %d]\n', M, N) ;

im_raw_resized = imresize(im_raw, sqrt(max_pixels / (M*N)) ) ;
[M, N] = size(im_raw_resized) ;
fprintf('image resized size : [%d, %d]\n', M, N) ;

im = single(im_raw_resized) ;

% frames to compute sift features
[p, q] = meshgrid(1:step:N, 1:step:M) ;
f = [p(:) q(:)]' ;
f(3, :) = 2 ;
f(4, :) = 0 ;

% sift at two bin_sizes
[f, d1] = vl_sift(im, 'frames', f, 'WindowSize', bin_size(1), 'verbose') ;
[f, d2] = vl_sift(im, 'frames', f, 'WindowSize', bin_size(2), 'verbose') ;

f = f(1:2, :) ;
d = vertcat(d1, d2) ;

fprintf('total number of features : %d\n', size(d, 2)) ;

end