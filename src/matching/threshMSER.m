function index = threshMSER(mserf,euclid_thresh)
% This function filters mser regions such that there are least overlapping
% regions. It takes mser regions array as input and returns the indices of
% filtered regions.

    %{
    clear all ;
    load('mser.mat') ;
    euclid_thresh = 4 ;
    mserf = mserf1 ;
    %}

    loc = mserf.Location ;
    sz = size(loc, 1) ;
    
    if sz <= 1
        index = 1:sz ;
        return
    end
    
    dis = pdist2(loc, loc) .^ 2 ;
    accept = ones(sz, 1) ;
    
    for i = 1 : sz
        if accept(i)==1
            xx = find(dis(i, :) < euclid_thresh) ;
            accept(xx) = 0 ;
            accept(i) = 1 ;
        end
    end
    index = find(accept==1) ;
end

%{
function index = threshMSER(mserf,thresh)
    loc1 = mserf.Location;
    loc2 = circshift(loc1,[1,0]);
    dist = sum((loc1-loc2).^2,2);
    index = find(dist>thresh);
end
%}
