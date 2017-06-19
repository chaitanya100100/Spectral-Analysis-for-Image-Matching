function I = normalize(I)
% this function linearly transform all values of input between 0 and 255.

mi = min(I(:)) ;
ma = max(I(:)) ;
I = 254 * (I - mi) / (ma - mi) ;
I = uint8(I) ;

end