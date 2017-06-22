function [f, d, M, N] = get_sift_features(im_path, max_pixels, step, bin_size)
% Resize image so that total number of feature points doesn't exceed much.
% There will be almost `max_pixels/(step*step)` feature points for small bin_size.
% Extract image SIFT features using vl_dsift.

bin_size = bin_size(1) ;

im = im2double(imread(im_path)) ;
if size(im, 3) == 3
    im = rgb2gray(im) ;
end

[M, N] = size(im) ;
im = imresize(im, sqrt(max_pixels / (M*N)) ) ;
im = im2single(im) ;
[M, N] = size(im) ;

% vl_dsift() does NOT compute a Gaussian scale space of the image I. 
% Instead, the image should be pre-smoothed at the desired scale level, e.b. by using the vl_imsmooth() function. 

[f, d] = vl_dsift(im, 'size', bin_size, 'step', step, 'verbose') ;

end