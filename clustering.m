%% Demo

% Parameters
params.p_data = [ pwd filesep 'dataset' filesep ] ;
params.libraries = [ 'matlab_bgl' ] ;
params.folders = [ 'clustering' ] ;
params.datasets = { 'PTC' , 'NCI109' , 'NCI1' , 'MUTAG' , 'ENZYMES' , 'DD' } ;

addpath( [ strrep(userpath,';','') filesep params.libraries ] ) ;
addpath( genpath([ pwd filesep params.folders ] ) ) ;

for dataset = 1:length(params.datasets)
    dataset_name = params.datasets{dataset} ;
    fprintf('Dataset = %s\n', dataset_name) ;
    
    % Load data
    data = load([params.p_data, dataset_name]) ;
    graphs = data.(dataset_name) ;
    labels = data.([ 'l' , lower(dataset_name) ]) ;
    
    % Randomly select one graph from each dataset to test the hierarchical
    % function
    graph_id = round(1 + (length(graphs)-1)*rand(1,1)) ;
    
    G = graphs(graph_id) ;
    lG = labels(graph_id) ;
    
    draw_graphs_from_adjmat( G.am ) ;
    
    % Construct hierarchy
    cluster = girvan_newman(G.am,floor(size(G.am,1)/2)) ;
    draw_graphs_clustering_from_adjmat(G.am, cluster);
    
    % We create an edge if we have any connection
    G1 = pyramid_construction(G, cluster, 0);
    draw_graphs_from_adjmat(G1.am);
    
    % connection ration 0.5
    G2 = pyramid_construction(G, cluster, 0.2);
    draw_graphs_from_adjmat(G2.am);
    
    cluster = grPartition(G.am,floor(size(G.am,1)/2)) ;
    draw_graphs_clustering_from_adjmat(G.am, cluster);
    
    % We create an edge if we have any connection
    G1 = pyramid_construction(G, cluster, 0);
    draw_graphs_from_adjmat(G1.am);
    
    % connection ration 0.5
    G2 = pyramid_construction(G, cluster, 0.2);
    draw_graphs_from_adjmat(G2.am);
    
    fprintf('Press any button...')
    pause;
    
end;