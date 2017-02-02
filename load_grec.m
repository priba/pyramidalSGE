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
        [graphs_train(i).v, graphs_train(i).nl, graphs_train(i).e,...
            graphs_train(i).el] = read_grec_gxl(fullfile(folder, fns_train{i})) ;
    end;

    for i = 1:nvalid
        [graphs_valid(i).v, graphs_valid(i).nl, graphs_valid(i).e,...
            graphs_valid(i).el] = read_grec_gxl(fullfile(folder, fns_valid{i})) ;
    end;

    for i = 1:ntest
        [graphs_test(i).v, graphs_test(i).nl, graphs_test(i).e,...
            graphs_test(i).el] = read_grec_gxl(fullfile(folder, fns_test{i})) ;
    end;

end