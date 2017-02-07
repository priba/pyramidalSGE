function [G1] = pyramid_construction(G , cluster, delta)
    G1 = G;
    % Create Nodes
    classes = unique(cluster);
    % Create Edges
    n_nodes = length(classes);
    G1.am = zeros(n_nodes, n_nodes) ;
    G1.nl.values = zeros(n_nodes, 1 );

    for i = 1:length(classes)
        idx_i = classes(i) == cluster ;
        
        % Discrete values
%        G1.nl.values(i) = round( mean(G.nl.values(idx_i,:),1) );        
        for j = i+1:length(classes)
            idx_j = classes(j) == cluster ;
            connections = G.am(idx_i, idx_j) ;
            ratio_connection = sum(connections(:)) / numel(connections) ;
            if ratio_connection > delta
                G1.am(i,j) = 1 ; G1.am(j,i) = 1 ;
            end ;
        end ;
    end ;
    G1.e = find(G1.am) ;
    G1.el.values = [] ;
end
