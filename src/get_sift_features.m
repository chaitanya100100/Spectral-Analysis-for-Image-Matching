function [f, d] = get_sift_features(im_path)
% Resize image so that total number of feature points doesn't exceed much.
% There will be almost `max_pixels/(step*step)` feature points for small bin_size.
% Extract image SIFT features using vlfeat.

max_pixels = 16000 ;
step = 6 ;
bin_size = 4 ;

im_raw = rgb2gray(imread('1.jpg')) ;
[M, N] = size(im_raw) ;
im_raw_resized = imresize(im_raw, sqrt(max_pixels / (M*N)) ) ;

im = single(im_raw_resized) ;

% vl_dsift() does NOT compute a Gaussian scale space of the image I. 
% Instead, the image should be pre-smoothed at the desired scale level, e.b. by using the vl_imsmooth() function. 

[f, d] = vl_dsift(im, 'size', bin_size, 'step', step, 'verbose') ;

end