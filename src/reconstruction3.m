function I = reconstruction(E, f, m, n)
% this function takes eigen vector values and put them at the feature
% points. Other values are interpolated linearly.

[X,Y] = meshgrid(1:n,1:m);
I = griddata(f(1,:)',f(2,:)', abs(E(:)), X, Y);

end