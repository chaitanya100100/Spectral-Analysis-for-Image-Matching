clear all ;

path1 = '../dataset/taj_day.jpg' ;
path2 = '../dataset/taj_night.jpg' ;

max_pixels = 16000 ;
step = 2 ;
bin_size = 4 ;
num = 5 ;

MinDiversity = 0.6 ;
MaxVariation = 0.3 ;
Delta = 15 ;


[II1, II2] = image_pair_eigen_spectram(path1, path2, max_pixels, step, bin_size, num) ;

I1 = II1{1} ;
I2 = II2{1} ;

[r1,f1] = mser(I1, MinDiversity, MaxVariation, Delta) ;
[r2,f2] = mser(I2, MinDiversity, MaxVariation, Delta) ;

figure ;
subplot(1, 2, 1) ;
imagesc(I1) ;
hold on ;
vl_plotframe(f1) ;

subplot(1, 2, 2) ;
imagesc(I2) ;
hold on ;
vl_plotframe(f2) ;