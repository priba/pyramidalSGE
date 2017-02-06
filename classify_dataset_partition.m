function [  ] = classify_dataset_partition( dataset_name, varargin )
%CLASSIFY_DATASET Summary of this function goes here
%   Detailed explanation goes here

    %% Parameters
        
    % Check hostname to know the file structure for adding the packages
    [~, name] = system('hostname');
    
    if ismember(upper(cellstr(name)), {'CVC206'} ) 
        % Dataset path
        params.p_data = '/home/adutta/Workspace/Datasets/STDGraphs' ;

        % External libreries (considering them to be on the userpath folder
        params.libraries_sub = {'matlab_bgl', 'libsvm/matlab',...
            'random_graphlet1'} ;
        
        % External libreries (considering them to be on the userpath
        % folder) No subfolders will be added
        params.libraries = { 'vlfeat/toolbox/vl_setup.m' } ;
        
        % Project folders
        params.folders = { 'clustering' } ;
        
        user_path = '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools' ;
        paths1 = pwd;
    elseif ismember(upper(cellstr(name)), {'CVC175'} )
        
        % Dataset path
        params.p_data = [ pwd filesep 'dataset' filesep ] ;

        % External libreries (considering them to be on the userpath
        % folder) All the subfolders will be added
        params.libraries_sub = { 'matlab_bgl', 'libsvm/matlab'};
        
        % External libreries (considering them to be on the userpath
        % folder) No subfolders will be added
        params.libraries = { 'vlfeat' } ;

        % Project folders
        params.folders = { 'clustering', 'random_graphlet1' } ;
        user_path = userpath;
        paths1 = pwd;
    else
        
        % Dataset path
        params.p_data = '/home/dag/adutta/Datasets/Graphs' ;

        % External libreries (considering them to be on the userpath folder
        params.libraries_sub = {'matlab_bgl', 'libsvm/matlab',...
            'random_graphlet1'} ;
        
        % External libreries (considering them to be on the userpath
        % folder) No subfolders will be added
        params.libraries = { 'vlfeat/toolbox/vl_setup.m' } ;
        
        % Project folders
        params.folders = { 'clustering' } ;
        
        user_path = '/home/dag/adutta/AdditionalTools' ;
        paths1 = '/home/dag/adutta/pyramidalSGE'
    end ;
    
    clear name;
    
    % Number of graphs depending on the number of edges
    params.T = [1 1 3 5 12 30 79 227 710 2322 8071]; 
    
    % Output folder
    params.out = [ pwd filesep 'out' filesep ];
    if ~exist( params.out , 'dir' )
        mkdir( params.out ) ;
    end ;
    
    % Output file
    params.headerSpec = 'Dataset: %s (%d graphs; %d classes; %d iterations)\n';
    params.sgeSpec = '\t*Stochastic Graphlet Embedding:\n\t\tEpsilon: %f\n\t\tDelta: %f\n';
    params.pyrSpec = '\t*Pyramidal:\n\t\tLevels: %d\n\t\tReduction: %f\n\t\tEdge Threshold: %f\n\t\tClustering function: %s\n';
    params.sepSpec = '----------------------------------------------------\n';

    %% Addpaths
    folders_paths = cellfun(@(x) strcat(paths1, filesep,x),params.folders, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath([folders_paths{:} ])

    user_path = strrep(user_path,';',''); % Windows
    user_path = strrep(user_path,':',''); % Linux
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries_sub, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath( [folders_paths{:} ])
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries, 'UniformOutput', false ) ;
    run( [folders_paths{:} ]) ;

    clear user_path folders_paths ;
    
    %% Default values
    [epsi, del, pyr_levels, pyr_reduction, edge_thresh, clustering_func, MAX2, node_label, nits, VERBOSE , task_id]...
        = input_parser( varargin ) ;
    rng(0);
        
    %% Loading
    if VERBOSE
        fprintf('Loading Dataset = %s\n', dataset_name) ;
    end
    [ graphs_train, clss_train, graphs_test, clss_test] ...
        =  load_database( params.p_data , dataset_name ) ;
    
    graphs = [graphs_train , graphs_test];
    clss = [clss_train; clss_test];
    %% Database information
    ntrain = size(graphs_train,2);
    ntest = size(graphs_test,2);
    ngraphs = size(graphs,2);
    classes = unique(clss);
    nclasses = size(classes,1);
    if VERBOSE
        fprintf('Dataset Information\n\tNumber of graph:%d\t(Train %d\tTest %d)\n\tNumber of classes: %d\n',...
            ngraphs,ntrain,ntest, nclasses) ;
    end
    
    %% Create histogram indices
  
    % Initialize storage
    for i = 1:pyr_levels
        M{i} = uint32(ceil(2*(params.T(1:MAX2(i))*log(2)+log(1/del))/epsi^2));
        global_var(i).hash_codes_uniq = cell(MAX2(i)-2,1);
        global_var(i).idx_image = cell(MAX2(i)-2,1);
        global_var(i).idx_bin = cell(MAX2(i)-2,1);
    end
    
    %% Iterate whole dataset
    for i = 1:ngraphs
        if VERBOSE
            fprintf('Graph: %d. ',i);
            tic;
        end ;
        
        pyr_graph = cell(pyr_levels, 1);
        pyr_graph{1} = graphs(i) ;
        
        % Construct hierarchy
        for j = 2:pyr_levels
            cluster = clustering_func(pyr_graph{j-1}.am,floor(size(pyr_graph{j-1}.am,1)/pyr_reduction)) ;
            pyr_graph{j} = pyramid_construction(pyr_graph{j-1}, cluster, edge_thresh);
        end;
        
        % Embedding
        for j = 1:pyr_levels
            if any(pyr_graph{j}.am(:))
                [ global_var(j) ] = graphlet_embedding(pyr_graph{j} , i , M{j} , global_var(j), MAX2(j) , node_label ) ;
            end
        end ;
        
        
        if VERBOSE
            toc
        end ;
    end;

    % Histogram dimensions
    dim_hists = cell(pyr_levels,1) ;
    for i = 1:pyr_levels
        dim_hists{i} = cellfun(@(x) size(x,1) ,global_var(i).hash_codes_uniq);
        clear global_var(i).hash_codes_uniq;
    end ;
    
    %% Compute histograms and kernels
    histograms = cell(pyr_levels,1);

    for j = 1:pyr_levels
        histograms{j} = cell(1,MAX2(j)-2);
        for i = 1:MAX2(j)-2
            histograms{j}{i} = sparse(global_var(j).idx_image{i},global_var(j).idx_bin{i},1,ngraphs,dim_hists{j}(i));
        end ;
    end ;
    
    % All possible combinations
    combinations = (1:MAX2(1)-2)';
    for j = 2:pyr_levels
        combinations = allcomb( { combinations, 1:MAX2(j)-2 } ) ;
    end ;
    
    maccs = zeros(size(combinations,1));
    
    for c = 1:size(combinations,1)
        
        comb = combinations(c,:);
        
        KM_train = zeros(ntrain,ntrain);
        KM_test = zeros(ntest,ntrain);

        % Concat histogram
        comb_hist = [];
        for i = 1:length(comb)
            comb_hist = [comb_hist , histograms{i}{comb(i)} ];
        end
        
        % Normalize hist
        comb_hist = comb_hist./repmat( sum(comb_hist,2)+eps ,1, size(comb_hist,2));
        
        X_train = comb_hist(1:ntrain,:);
        X_test = comb_hist(ntrain+(1:ntest),:);

        KM_train(:,:) = vl_alldist2(X_train',X_train','KL1');
        KM_test(:,:) = vl_alldist2(X_test',X_train','KL1');


        %% Evaluate
        % Evaluate nits times to get the accuracy mean and standard deviation
        train_classes = clss(1:ntrain);
        test_classes = clss(ntrain+(1:ntest));


        % Training and testing individual kernels

        K_train = [(1:ntrain)' KM_train];
        K_test = [(1:ntest)' KM_test];

        cs = 5:5:100;
        best_cv = 0;

        for j = 1:length(cs)

            options = sprintf('-s 0 -t 4 -v %d -c %f -b 1 -g 0.07 -h 0 -q',...
                    nits,cs(j));
            model_libsvm = svmtrain(train_classes,K_train,options);

            if(model_libsvm>best_cv)
                best_cv = model_libsvm;
                best_c = cs(j);
            end;

        end;

        options = sprintf('-s 0 -t 4 -c %f -b 1 -g 0.07 -h 0 -q',...
            best_c);

        model_libsvm = svmtrain(train_classes,K_train,options);

        [~,acc,~] = svmpredict(test_classes,K_test,model_libsvm,'-b 1');

        % Mean
        maccs(c) = acc(1);
        
    end
    clear global_var;
    
    % Save results
    if task_id<0
        fileID = fopen([params.out dataset_name '_' node_label '.txt'],'a') ;
    else
        fileID = fopen([params.out dataset_name '_' node_label '_' num2str(task_id) '.txt'],'a') ;
    end
    fprintf(fileID,params.headerSpec, dataset_name, ngraphs, nclasses, nits) ;
    fprintf(fileID,params.sgeSpec, epsi , del) ;
    fprintf(fileID,params.pyrSpec, pyr_levels , pyr_reduction , edge_thresh , func2str(clustering_func)) ;
    for i = 1:size(combinations,1)
        fprintf(fileID, 't = ');
        if VERBOSE
            fprintf('t = ');
        end ;
        for j = 1:size(combinations,2)
            fprintf(fileID, '%d\t', combinations(i,j)+2) ;
            if VERBOSE
            	fprintf('%d\t', combinations(i,j)+2) ;
            end ;
        end ;
        fprintf(fileID, '%.2f\n', maccs(i));
        if VERBOSE
            fprintf('%.2f\n', maccs(i));
        end ;
    end;
    fprintf(fileID,params.sepSpec) ;
    fclose(fileID);
    
    
    %% Rmpaths
    folders_paths = cellfun(@(x) strcat(pwd, filesep,x),params.folders, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    rmpath([folders_paths{:} ])
    
    user_path = strrep(userpath,';',''); % Windows
    user_path = strrep(user_path,':',''); % Linux
    

    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries_sub, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    rmpath( [folders_paths{:} ])
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries, 'UniformOutput', false ) ;
    rmpath( [folders_paths{:} ])
    
    clear user_path folders_paths ;
    
end

%% Loads the database in the correct format
function [ graphs_train, clss_train, graphs_test, clss_test] = ...
    load_database( p_data , dataset_name )
    switch dataset_name
        case 'GREC'
            [ graphs_train, clss_train, graphs_valid, clss_valid,...
                graphs_test, clss_test] = load_grec([ p_data filesep dataset_name filesep 'data']);
        case 'GWHistoGraph'
            [ graphs_train, clss_train, graphs_valid, clss_valid,...
                graphs_test, clss_test] = load_gw(p_data);            
        otherwise
            error('pyramidalSGE:incorrectDataset',...
                'Error. \nDatabase %s not accepted.', dataset_name)
    end;
    graphs_train = [graphs_train , graphs_valid];
    clss_train = [clss_train ; clss_valid];
end
