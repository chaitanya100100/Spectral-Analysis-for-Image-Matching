clear all ;

path1 = '1.jpg' ;
path2 = '2.jpg' ;

fprintf('dense sift features for %s\n', path1) ;
[f1, d1] = get_sift_features(path1) ;

fprintf('dense sift features for %s\n', path2) ;
[f2, d2] = get_sift_features(path2) ;

