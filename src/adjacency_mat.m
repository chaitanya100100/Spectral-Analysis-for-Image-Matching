function A = adjacency_mat(d1, d2)
% this function creates adjacency matrix withing the vertices of each image 
% and across the vertices of two images and thus making complete adjacency matrix

d1 = double(d1') ;
d2 = double(d2') ;
Adj1 = squareform(pdist(d1, 'cosine')) ;
Adj2 = squareform(pdist(d2, 'cosine')) ;
Adjcross = pdist2(d1, d2, 'cosine') ;

A_upper = horzcat(Adj1, Adjcross) ;
A_lower = horzcat(Adjcross', Adj2) ;
A = vertcat(A_upper, A_lower) ;
A = A.*A ;
A = exp(-A) ;
end