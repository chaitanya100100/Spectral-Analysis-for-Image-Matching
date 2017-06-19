% This script first loads a mat file which has eigen function computed for an
% image pair and then extract mser features from each of them. Then it
% matches extracted features and show matched features.

clear all;
addpath('./modified') ;

fold = 'madrid' ;
load(strcat('../results/', fold, '/', fold,'.mat')) ;

num = size(IS1, 2) ;
im1 = imread(strcat('../dataset/', fold, '/01.png')) ;
im1 = imresize(im1, size(IS1{1})) ;

im2 = imread(strcat('../dataset/', fold, '/02.png')) ;
im2 = imresize(im2, size(IS2{1})) ;

for i = 1 : num

I1 = IS1{i} ;
I2 = IS2{i} ;

I1 = normalize(I1) ;
I2 = normalize(I2) ;

points1 = detectMSERFeatures(I1, 'MaxAreaVariation', 0.5, 'ThresholdDelta', 2);
points2 = detectMSERFeatures(I2, 'MaxAreaVariation', 0.5, 'ThresholdDelta', 2); 


[f1, vpts1] = extractFeatures(I1, points1, 'SURFSize', 64) ;
[f2, vpts2] = extractFeatures(I2, points2, 'SURFSize', 64) ;

indexPairs = matchFeatures(f1, f2, 'Unique', true,'MaxRatio',0.7,'MatchThreshold',20) ;
matched_pts1 = vpts1(indexPairs(:, 1)) ;
matched_pts2 = vpts2(indexPairs(:, 2)) ;


matched_pts1 = matched_pts1.Location ;
matched_pts2 = matched_pts2.Location ;

figure;
stackedImage = cat(2, I1, I2); % Places the two images side by side
imagesc(stackedImage); colormap jet;
width = size(im1, 2);
hold on;
numPoints = size(matched_pts1, 1); % points2 must have same # of points
% Note, we must offset by the width of the image
for i = 1 : numPoints
    plot(matched_pts1(i, 1), matched_pts1(i, 2), 'y+', matched_pts2(i, 1) + width, ...
         matched_pts2(i, 2), 'y+');
    line([matched_pts1(i, 1) matched_pts2(i, 1) + width], [matched_pts1(i, 2) matched_pts2(i, 2)], ...
         'Color', 'black');
end

end

%{
figure; showMatchedFeatures(im1,im2,matched_pts1,matched_pts2,'montage');
legend('matched points 1','matched points 2');
%}


%{
figure ;
imagesc(I1) ; colormap jet ;
hold on ; plot(points1) ;

figure ;
imagesc(I2) ; colormap jet ;
hold on ; plot(points2) ;
%}
