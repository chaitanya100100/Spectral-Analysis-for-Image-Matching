function js = jspec(f, I) 
    
    mul = 4 ;
    [M, N] = size(I) ;
    js = zeros(size(f, 2), 128) ;
    
    for i = 1 : size(f, 2)
        x = f(:, i) ;
        x(3) = sqrt(x(3))*mul ;
        x(5) = sqrt(x(5))*mul ;
        x = uint16(x) ;
        part = I(max(1, x(2) - x(5)) : min(M, x(2) + x(5)), max(1, x(1) - x(3)) : min(N, x(1) + x(3))) ;
        p = max(size(part)) ;
        part = imresize(part, [p, p]) ;
        [~, js(i, :)] = vl_dsift(im2single(part), 'size', floor(p/4), 'step', p) ; 
    end
    js = horzcat(js, f') ;
end