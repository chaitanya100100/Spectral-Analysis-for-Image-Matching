% This script loads mat file which has eigen functions of an image pair. It
% then extract mser features, computes jspec features and then runs knn
% search with some constraints for matching i.e. bidirectional matching.

clear all ;
load('montreal.mat') ;

MinDiversity = 0.9 ;
MaxVariation = 0.4 ;
Delta = 10 ;
num = size(IS1, 2) ;


% extract mser features
for i = 2:num
    I1 = IS1{i} ;
    I2 = IS2{i} ;
    [r1,f1,I1] = mser(I1, MinDiversity, MaxVariation, Delta) ;
    [r2,f2,I2] = mser(I2, MinDiversity, MaxVariation, Delta) ;
    
    js1 = jspec(f1, I1) ;
    js2 = jspec(f2, I2) ;
    
    [nn12, dist12] = knnsearch(js2, js1, 'K', 1) ;
    [nn21, dist21] = knnsearch(js1, js2, 'K', 1) ;
    
    k = 0 ;
    for j = 1 : size(nn12)
        if nn21(nn12(j)) == j
            k = k + 1 ;
            nf1(:, k) = f1(:, j) ;
            nf2(:, k) = f2(:, nn12(j)) ;
        end
    end
    
    
    figure ;
    subplot(1, 2, 1) ;
    imagesc(I1) ; colormap jet ;
    hold on ;
    vl_plotframe(nf1) ;

    subplot(1, 2, 2) ;
    imagesc(I2) ; colormap jet ;
    hold on ;
    vl_plotframe(nf2) ;
    
end
