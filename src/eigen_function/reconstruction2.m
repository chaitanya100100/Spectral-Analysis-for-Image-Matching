function I = reconstruction2(E, f, m, n)
% this function takes eigen vector values and put them at the feature
% points. Other values are kept zero i.e. no interpolation.

%{
c = unique(f(1, :)) ;
r = unique(f(2, :)) ;
[C, R] = meshgrid(c, r) ;
c = size(c, 2) ;
r = size(r, 2) ;
%}

I = zeros(m, n) ;
I(sub2ind(size(I), f(2, :), f(1, :))) = E(:) ;

end