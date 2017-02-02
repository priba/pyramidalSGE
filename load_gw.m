function [ graphs_train, clss_train, graphs_valid, clss_valid, ...
    graphs_test, clss_test] = load_gw(folder, subfolder)
    %% Load GWHistoGraph database
    fold = '/Data/Word_Graphs/01_Skew' ;

    [fns_train,clss_train] = read_gw_txt(fullfile(folder, 'Set','Train.txt')) ;
    [fns_valid,clss_valid] = read_gw_txt(fullfile(folder, 'Set','Valid.txt')) ;
    [fns_test,clss_test] = read_gw_txt(fullfile(folder, 'Set', 'Test.txt')) ;

    ntrain = length(fns_train) ; 
    nvalid = length(fns_valid) ; 
    ntest = length(fns_test) ; 

    for i = 1:ntrain
        [graphs_train(i).v, graphs_train(i).nl, graphs_train(i).e,...
            graphs_train(i).el] = read_gw_gxl(fullfile(folder, fold,...
            subfolder, fns_train{i})) ;
    end;

    for i = 1:nvalid
        [graphs_valid(i).v, graphs_valid(i).nl, graphs_valid(i).e,...
            graphs_valid(i).el] = read_gw_gxl(fullfile(folder, fold,...
            subfolder, fns_valid{i})) ;
    end;

    for i = 1:ntest
        [graphs_test(i).v, graphs_test(i).nl, graphs_test(i).e,...
            graphs_test(i).el] = read_gw_gxl(fullfile(folder, fold,...
            subfolder, fns_test{i})) ;
    end;

end