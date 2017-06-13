% run('/home/chaitanya/CVIT/2017-05-23_summer/vlfeat/vlfeat-0.9.20/toolbox/vl_setup.m')
clear all ;

path1 = '../dataset/01.png' ;
path2 = '../dataset/02.png' ;

max_pixels = 32000 ;
step = 5 ;
bin_size = [6 10] ;
num = 5 ;

MinDiversity = 0.5 ;
MaxVariation = 0.3 ;
Delta = 10 ;


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

IS1 = {} ;
IS2 = {} ;
% reconstruction of spectral images
for i = 1 : num
    E1 = V(1:nf1, i) ;
    E2 = V(nf1+1:nf, i) ;

    IS1{i} = reconstruction(E1, f1, m1, n1) ;
    IS2{i} = reconstruction(E2, f2, m2, n2) ;
end

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
%{
% extract mser features
[r1,f1,I1] = mser(I1, MinDiversity, MaxVariation, Delta) ;
[r2,f2,I2] = mser(I2, MinDiversity, MaxVariation, Delta) ;

fr1 = f1(1:2, :) ;
fr1(3, :) = 1 ;
fr1(4, :) = 0 ;
[sf1, sd1] = vl_sift(single(I1), 'frames', fr1) ;


fr2 = f2(1:2, :) ;
fr2(3, :) = 1 ;
fr2(4, :) = 0 ;
[sf2, sd2] = vl_sift(single(I2), 'frames', fr2) ;

[nn12, dist12] = knnsearch(sd2', sd1', 'K', 2) ;
[nn21, dist21] = knnsearch(sd1', sd2', 'K', 2) ;


figure ;
subplot(1, 2, 1) ;
imagesc(I1) ;
hold on ;
vl_plotframe(f1) ;

subplot(1, 2, 2) ;
imagesc(I2) ;
hold on ;
vl_plotframe(f2) ;
%}
