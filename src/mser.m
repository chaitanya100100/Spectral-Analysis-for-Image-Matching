function [r, f, I] = mser(I, MinDiversity, MaxVariation, Delta)

mi = min(I(:)) ;
ma = max(I(:)) ;
I = 255 * (I - mi) / (ma - mi) ;
I = uint8(I) ;

[r,f] = vl_mser(I,'MinDiversity',MinDiversity,'MaxVariation',MaxVariation,'Delta',Delta) ;
f = vl_ertr(f) ;

end