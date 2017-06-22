function [r, f, I] = mser(I, MinDiversity, MaxVariation, Delta)
% This function returns mser regions calculated by vl_mser function after
% normalizing the image.

mi = min(I(:)) ;
ma = max(I(:)) ;
I = 254 * (I - mi) / (ma - mi) ;
I = uint8(I) ;

[r,f] = vl_mser(I,'MinDiversity',MinDiversity,'MaxVariation',MaxVariation,'Delta',Delta) ;
f = vl_ertr(f) ;

end