clear all ;

fold = 'stargarder' ;

addpath('../modified') ;

load(strcat('../../results/', fold, '/', fold,'.mat')) ;
num = size(IS1, 2) ;

for j = 1:num

I1 = IS1{j} ;
I2 = IS2{j} ;
%I2 = imresize(I2, size(I1)) ;

I1 = normalize(I1) ;
I2 = normalize(I2) ;

mserf1 = detectMSERFeatures(I1, 'MaxAreaVariation', 0.5, 'ThresholdDelta', 1) ;
mserf2 = detectMSERFeatures(I2, 'MaxAreaVariation', 0.5, 'ThresholdDelta', 1) ; 

mserf1 = mserf1(threshMSER(mserf1, 5)) ;
mserf2 = mserf2(threshMSER(mserf2, 5)) ;

loc1 = mserf1.Location ;
loc2 = mserf2.Location ;

[surff1, vpts1] = extractFeatures(I1, loc1, 'Method', 'SURF') ;
[surff2, vpts2] = extractFeatures(I2, loc2, 'Method', 'SURF') ;

indexPairs = match_features(mserf1, mserf2, surff1, surff2) ;

matched_pts1 = vpts1(indexPairs(:, 1)) ;
matched_pts2 = vpts2(indexPairs(:, 2)) ;

matched_pts1 = matched_pts1.Location ;
matched_pts2 = matched_pts2.Location ;

fx = figure; axis off ;
stackedImage = cat(2, I1, I2); % Places the two images side by side
imagesc(stackedImage); colormap jet;
width = size(I1, 2);
hold on;
numPoints = size(matched_pts1, 1); % points2 must have same # of points
% Note, we must offset by the width of the image
for i = 1 : numPoints
    plot(matched_pts1(i, 1), matched_pts1(i, 2), 'y+', matched_pts2(i, 1) + width, ...
         matched_pts2(i, 2), 'y+');
    line([matched_pts1(i, 1) matched_pts2(i, 1) + width], [matched_pts1(i, 2) matched_pts2(i, 2)], ...
         'Color', 'black');
end
axis off ;
saveas(fx, strcat(num2str(j), '.jpg')) ;

end