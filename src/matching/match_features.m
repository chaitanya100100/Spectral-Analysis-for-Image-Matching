function features = match_features(mserf1,mserf2,surff1,surff2)
    thresh_euclid=3;
    mserthresh = 0.6;
    surfthresh = 0.6;
    
    sigma = 1;
    surfdist = pdist2(surff1,surff2).^2;
    mserdist = pdist2(mserf1.Location,mserf2.Location).^2;
    [val,idx] = sort(surfdist,2);
    accept = zeros(size(idx,1),1);
    for i=1:size(idx)
        if(mserdist(i,idx(i,1))<mserthresh*mserdist(i,idx(i,2)) & ...
               surfdist(i,idx(i,1))<surfthresh*surfdist(i,idx(i,2)))
            accept(i)=1;
        end
    end
   
    ind = find(accept==1);
    features(:,1)=ind;
    features(:,2)=idx(ind,1);
end
