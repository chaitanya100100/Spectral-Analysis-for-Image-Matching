% This script first loads a mat file which has eigen function computed for an
% image pair and then extract mser features from each of them. Finally it
% displays detected features.

clear all;

fold = 'madrid' ;
load(strcat('../../results/', fold, '/', fold,'.mat')) ;

MinDiversity = 0.9 ;
MaxVariation = 0.5 ;
Delta = 8 ;
num = size(IS1, 2) ;

x = figure ;
%subplot(1, 2, 1) ;
I1 = imread(strcat('../dataset/', fold, '/01.png')) ;
I1 = imresize(I1, size(IS1{1})) ;
imshow(I1) ;
%saveas(x, strcat('../results/', fold, '/r1.png') ) ;

x = figure ;
%subplot(1, 2, 2) ;
I2 = imread(strcat('../dataset/', fold, '/02.png')) ;
I2 = imresize(I2, size(IS2{1})) ;
imshow(I2) ;
%saveas(x, strcat('../results/', fold, '/r2.png') ) ;

% extract mser features
for i = 2:num
    I1 = IS1{i} ;
    I2 = IS2{i} ;
    [r1,f1,I1] = mser(I1, MinDiversity, MaxVariation, Delta) ;
    [r2,f2,I2] = mser(I2, MinDiversity, MaxVariation, Delta) ;
    
    x = figure ;
    %subplot(1, 2, 1) ;
    imagesc(I1) ; colormap jet ; axis equal off ;
    hold on ;
    vl_plotframe(f1) ; axis tight ;
    %saveas(x, strcat('../results/', fold, '/e', num2str(i),'_1.png') ) ;
    
    x = figure ;
    %subplot(1, 2, 2) ;
    imagesc(I2) ; colormap jet ; axis equal off ;
    hold on ;
    vl_plotframe(f2) ; axis tight ;
    %saveas(x, strcat('../results/', fold, '/e', num2str(i),'_2.png') ) ;
    %print(x, '-dpng', strcat('../results/', fold, '/e', num2str(i),'.png'), '-bestfit' ) ;
    
end
