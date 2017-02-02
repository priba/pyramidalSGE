%% Run classification script
function [] = run_sge(dataset_name)

clc;

rng(0); 

addpath(genpath('./clustering'));

% Dataset
% dataset_name = 'GREC';

% Information
VERBOSE = 1 ;

% Embedding
eps_i = [ 0.1 , 0.05 ] ;
del_i = [ 0.1 , 0.05 ] ;
max2 = [7, 5, 5];
node_label = 'label';
% Pyramid
pyr_levels = [ 2 , 1 , 2  ] ;
pyr_reduction = 2 ;
edge_tresh = 0 ;
clustering_func = @girvan_newman ;

% Standard error
nits = 10 ;

for eps = eps_i
    for del = del_i
        for pyr_level = pyr_levels
            if ismember(lower(dataset_name), {'grec','gwhistograph'})
                classify_dataset_partition(dataset_name, 'VERBOSE', VERBOSE, 'epsilon', eps, 'delta', del, ...
                    'pyr_levels', pyr_level,'pyr_reduction', pyr_reduction, 'edge_tresh', edge_tresh, ...
                    'max2', max2(1:pyr_level), 'label', node_label, 'clustering_func' , clustering_func);
            else
                classify_dataset_kfold(dataset_name, 'VERBOSE', VERBOSE, 'epsilon', eps, 'delta', del, ...
                    'pyr_levels', pyr_level,'pyr_reduction', pyr_reduction, 'edge_tresh', edge_tresh, ...
                    'max2', max2(1:pyr_level), 'label', node_label, 'clustering_func' , clustering_func);
            end
        end ;
    end ;
end ;