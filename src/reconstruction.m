function I = reconstruction(E, f, m, n)
% this function takes eigen vector values and put them at the feature
% points. Other values are interpolated linearly.

%{
c = unique(f(1, :)) ;
r = unique(f(2, :)) ;
[C, R] = meshgrid(c, r) ;
c = size(c, 2) ;
r = size(r, 2) ;
%}


c = size(unique(f(1, :)), 2) ;
r = size(unique(f(2, :)), 2) ;
C = reshape(f(1, :), [r c]) ;
R = reshape(f(2, :), [r c]) ;


ER = reshape(E, [r c]) ;
[NC, NR] = meshgrid(1:n, 1:m) ;

I = interp2(C, R, abs(ER), NC, NR, 'linear', 0) ;
end