function [ graphs_train, clss_train, graphs_valid, clss_valid,...
    graphs_test, clss_test] = load_grec(folder)
    %% Load GREC database

    [fns_train,clss_train] = read_grec_cxl(fullfile(folder,'train.cxl')) ;
    [fns_valid,clss_valid] = read_grec_cxl(fullfile(folder,'valid.cxl')) ;
    [fns_test,clss_test] = read_grec_cxl(fullfile(folder,'test.cxl')) ;

    ntrain = length(fns_train) ; 
    nvalid = length(fns_valid) ; 
    ntest = length(fns_test) ; 

    for i = 1:ntrain
        [graphs_train(i).v, graphs_train(i).nl.values, graphs_train(i).e,...
            graphs_train(i).el.values] = read_grec_gxl(fullfile(folder, fns_train{i})) ;
        graphs_train(i).am = zeros(size(graphs_train(i).v,1),size(graphs_train(i).v,1)) ;
        ind = sub2ind(size(graphs_train(i).am),graphs_train(i).e(:,1),graphs_train(i).e(:,2)) ;
        graphs_train(i).am(ind) = 1 ;
%         graphs_train(i).am = double(graphs_train(i).am |graphs_train(i).am');
    end;

    for i = 1:nvalid
        [graphs_valid(i).v, graphs_valid(i).nl.values, graphs_valid(i).e,...
            graphs_valid(i).el.values] = read_grec_gxl(fullfile(folder, fns_valid{i})) ;
        graphs_valid(i).am = zeros(size(graphs_valid(i).v,1),size(graphs_valid(i).v,1)) ;
        ind = sub2ind(size(graphs_valid(i).am),graphs_valid(i).e(:,1),graphs_valid(i).e(:,2)) ;
        graphs_valid(i).am(ind) = 1 ;
%         graphs_valid(i).am = double(graphs_valid(i).am |graphs_valid(i).am');
    end;

    for i = 1:ntest
        [graphs_test(i).v, graphs_test(i).nl.values, graphs_test(i).e,...
            graphs_test(i).el.values] = read_grec_gxl(fullfile(folder, fns_test{i})) ;
        graphs_test(i).am = zeros(size(graphs_test(i).v,1),size(graphs_test(i).v,1)) ;
        ind = sub2ind(size(graphs_test(i).am),graphs_test(i).e(:,1),graphs_test(i).e(:,2)) ;
        graphs_test(i).am(ind) = 1 ;
%         graphs_test(i).am = double(graphs_test(i).am |graphs_test(i).am');
    end;

end