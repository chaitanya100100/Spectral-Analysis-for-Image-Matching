%run('/home/chaitanya/CVIT/2017-05-23_summer/vlfeat/vlfeat-0.9.20/toolbox/vl_setup.m')

function [IS1, IS2] = main(fold)
% this function generates eigen functions for image pair and save it to a
% mat file. To visualize it, run 'see_matching.m' script.

%fold = 'montreal' ;

path1 = strcat('../dataset/', fold, '/01.png') ;
path2 = strcat('../dataset/', fold, '/02.png') ;

max_pixels = 96000 ;
step = 4 ;
bin_size = [10 6] ;
num = 6 ;

MinDiversity = 0.9 ;
MaxVariation = 0.4 ;
Delta = 10 ;


% extract dense sift features from both images
fprintf('dense sift features for %s\n', path1) ;
[f1, d1, m1, n1] = get_sift_features_128D(path1, max_pixels, step, bin_size) ;
fprintf('dense sift features for %s\n', path2) ;
[f2, d2, m2, n2] = get_sift_features_128D(path2, max_pixels, step, bin_size) ;

% number of feature points of image1 and image2 is nf1 and nf2 respectively
nf1 = size(d1, 2) ;
nf2 = size(d2, 2) ;
nf = nf1 + nf2 ;

% joint graph spectram
A = adjacency_mat(d1, d2) ;
[V, E] = graph_eigen_spectram(A, num) ;

IS1 = {} ;
IS2 = {} ;
% reconstruction of spectral images
for i = 1 : num
    E1 = V(1:nf1, i) ;
    E2 = V(nf1+1:nf, i) ;

    IS1{i} = reconstruction(E1, f1, m1, n1) ;
    IS2{i} = reconstruction(E2, f2, m2, n2) ;
end

save(strcat(fold, '.mat'), 'IS1', 'IS2') ;

%{
for i = 2 : num
    I1 = IS1{i} ;
    I2 = IS2{i} ;

    figure ;
    subplot(1, 2, 1) ;
    imagesc(I1) ;
    colormap jet ;
    subplot(1, 2, 2) ;
    imagesc(I2) ;
    colormap jet ;
end
%}

%{
% extract mser features
for i = 2:num
    I1 = IS1{i} ;
    I2 = IS2{i} ;
    [r1,f1,I1] = mser(I1, MinDiversity, MaxVariation, Delta) ;
    [r2,f2,I2] = mser(I2, MinDiversity, MaxVariation, Delta) ;

    figure ;
    subplot(1, 2, 1) ;
    imagesc(I1) ; colormap jet ;
    hold on ;
    vl_plotframe(f1) ;

    subplot(1, 2, 2) ;
    imagesc(I2) ; colormap jet ;
    hold on ;
    vl_plotframe(f2) ;
end
%}

%{
fr1 = f1(1:2, :) ;
fr1(3, :) = 10 ;
fr1(4, :) = 0 ;
[sf1, sd1] = vl_sift(single(I1), 'frames', fr1) ;


fr2 = f2(1:2, :) ;
fr2(3, :) = 10 ;
fr2(4, :) = 0 ;
[sf2, sd2] = vl_sift(single(I2), 'frames', fr2) ;

[nn12, dist12] = knnsearch(sd2', sd1', 'K', 2) ;
[nn21, dist21] = knnsearch(sd1', sd2', 'K', 2) ;
%}
end