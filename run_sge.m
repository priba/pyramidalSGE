addpath(genpath('./clustering'));

VERBOSE = 0 ;
eps_i = [ 0.1 , 0.05 ] ;
del_i = [ 0.1 , 0.05 ] ;
pyr_levels = [ 1 , 2 ] ;
pyr_reduction = 2 ;
edge_tresh = 0 ;
clustering_func = @girvan_newman ;
nits = 10 ;

for eps = eps_i
    for del = del_i
        for pyr_level = pyr_levels
            classify_dataset('MUTAG', 'VERBOSE',VERBOSE, 'epsilon', eps, 'delta', del, ...
                'pyr_levels',pyr_level,'pyr_reduction',pyr_reduction, 'edge_tresh',edge_tresh, ...
                'clustering_func' , clustering_func);
        end ;
    end ;
end ;