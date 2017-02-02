function [  ] = classify_dataset( dataset_name, varargin )
%CLASSIFY_DATASET Summary of this function goes here
%   Detailed explanation goes here

    %% Parameters
        
    % Check hostname to know the file structure for adding the packages
    [~, name] = system('hostname');
    
    if ismember(upper(cellstr(name)), {'CVC206'} ) 
        % Dataset path
        params.p_data = '/home/adutta/Workspace/Datasets/STDGraphs' ;

        % External libreries (considering them to be on the userpath folder
        params.libraries = {'matlab_bgl', 'libsvm/matlab', 'vlfeat',...
            'random_graphlet1'} ;

             % Project folders
        params.folders = { 'clustering' } ;
        
        user_path = '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools' ;
    else
        % Dataset path
        params.p_data = [ pwd filesep 'dataset' filesep ] ;

        % External libreries (considering them to be on the userpath folder
        params.libraries = { 'matlab_bgl', 'libsvm/matlab', 'vlfeat' } ;

        % Project folders
        params.folders = { 'clustering', 'random_graphlet1' } ;
        user_path = userpath;
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
    folders_paths = cellfun(@(x) strcat(pwd, filesep,x),params.folders, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath([folders_paths{:} ])
    
    user_path = strrep(user_path,';',''); % Windows
    user_path = strrep(user_path,':',''); % Linux
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath( [folders_paths{:} ])

    clear user_path folders_paths ;
    
    %% Default values
    [epsi, del, pyr_levels, pyr_reduction, edge_thresh, clustering_func, MAX2, node_label, nits, VERBOSE]...
        = input_parser( varargin ) ;
    rng(0);
        
    %% Loading
    if VERBOSE
        fprintf('Loading Dataset = %s\n', dataset_name) ;
    end
    [ graphs , clss ] =  load_database( params.p_data , dataset_name ) ;

    %% Database information
    ngraphs = size(graphs,2);
    classes = unique(clss);
    nclasses = size(classes,1);
    if VERBOSE
        fprintf('Dataset Information\n\tNumber of graph: %d\n\tNumber of classes: %d\n',...
            ngraphs, nclasses) ;
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
        combinations = allcomb(combinations, 1:MAX2(j)-2) ;
    end ;
    
    maccs = zeros(size(combinations,1));
    mstds = zeros(size(combinations,1));
    for c = 1:size(combinations,1)
        
        comb = combinations(c,:);
        
        KM_train = zeros(ngraphs,ngraphs);
        KM_test = zeros(ngraphs,ngraphs);

        % Concat histogram
        comb_hist = [];
        for i = 1:length(comb)
            comb_hist = [comb_hist , histograms{i}{comb(i)} ];
        end
        
        % Normalize hist
        comb_hist = comb_hist./repmat( sum(comb_hist,2)+eps ,1, size(comb_hist,2));
        
        X_train = comb_hist;
        X_test = comb_hist;

        KM_train(:,:) = vl_alldist2(X_train',X_train','KL1');
        KM_test(:,:) = vl_alldist2(X_test',X_train','KL1');


        %% Evaluate
        % Evaluate nits times to get the accuracy mean and standard deviation
        accs = zeros(nits,1);
        for it = 1:nits
            train_set = [];

            for i = classes'
                idx = find(clss == i);
                lngth_idx = length(idx);
                p90 = round(lngth_idx*0.9);
                train_set = [train_set;randsample(idx,p90)];

            end;

            train_set = sort(train_set);
            test_set = setdiff(1:ngraphs,train_set')';

            train_classes = clss(train_set);
            test_classes = clss(test_set);

            ntrain_set = size(train_set,1);
            ntest_set = size(test_set,1);

            % Training and testing individual kernels

            K_train = [(1:ntrain_set)' KM_train(train_set,train_set)];
            K_test = [(1:ntest_set)' KM_test(test_set,train_set)];

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
            accs(it) = acc(1);
        end ;

        % Mean and standard deviation
        maccs(c) = mean(accs);
        mstds(c) = std(accs)./sqrt(nits);
    end
    clear global_var;
    
    % Save results
    fileID = fopen([params.out dataset_name '_' node_label '.txt'],'a') ;
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
        fprintf(fileID, '%.2f \\pm %.2f \n', maccs(i),mstds(i));
        if VERBOSE
            fprintf('%.2f \\pm %.2f \n', maccs(i),mstds(i));
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
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    rmpath( [folders_paths{:} ])

    clear user_path folders_paths ;
    
end

%% Loads the database in the correct format
function [ graphs , clss ] = load_database( p_data , dataset_name )
    switch dataset_name
        case 'MUTAG'
            load(fullfile(p_data,'MUTAG.mat'));
            lmutag(lmutag == -1) = 2;clss = lmutag;
            graphs = MUTAG;
        case 'PTC'
            load(fullfile(p_data,'PTC.mat'));
            lptc(lptc == -1) = 2;clss = lptc;
            graphs = PTC;
        case 'ENZYMES'
            load(fullfile(p_data,'ENZYMES.mat'));
            clss = lenzymes;
            graphs = ENZYMES;
        case 'DD'
            load(fullfile(p_data,'DD.mat'));
            clss = ldd;
            graphs = DD;
        case 'NCI1'
            load(fullfile(p_data,'NCI1.mat'));
            clss = lnci1;
            graphs = NCI1;
        case 'NCI109'
            load(fullfile(p_data,'NCI109.mat'));
            clss = lnci109;
            graphs = NCI109;
        case 'GREC'
            [ graphs_train, clss_train, graphs_valid, clss_valid,...
                graphs_test, clss_test] = load_grec([ p_data filesep 'data']);
        case 'GWHistoGraph'
            [ graphs_train, clss_train, graphs_valid, clss_valid,...
                graphs_test, clss_test] = load_gw(p_data, subfolder);            
        otherwise
            error('pyramidalSGE:incorrectDataset',...
                'Error. \nDatabase %s not accepted.', dataset_name)
    end;
end

%% Load GREC database
function [ graphs_train, clss_train, graphs_valid, clss_valid,...
    graphs_test, clss_test] = load_grec(folder)

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

%% Load GWHistoGraph database
function [ graphs_train, clss_train, graphs_valid, clss_valid, ...
    graphs_test, clss_test] = load_gw(folder, subfolder)

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

%% Default values
function [eps, del, pyr_levels, pyr_reduction, edge_tresh, clustering_func,...
    max2, node_label, nits, VERBOSE] = input_parser( input )
    VERBOSE = 0 ;
    eps = 0.1 ;
    del = 0.1 ;
    pyr_levels = 1 ;
    pyr_reduction = 2 ;
    edge_tresh = 0 ;
    clustering_func = @girvan_newman ;
    
    node_label = 'unlabel' ;
    
    nits = 10 ;
    
    % Parse optional input parameters
    v = 1;
    while v < numel(input)
        switch lower(input{v})
        case 'verbose'
            v = v+1;
            VERBOSE = input{v};
        case 'epsilon'
            v = v+1;
            eps = input{v};
        case 'delta'
            v = v+1;
            del = input{v};
        case 'pyr_levels'
            v = v+1;
            pyr_levels = input{v};
        case 'pyr_reduction'
            v = v+1;
            pyr_reduction = input{v};
            assert(pyr_reduction>=1) ;
        case 'edge_tresh'
            v = v+1;
            edge_tresh = input{v};
            assert(edge_tresh<=1 && edge_tresh>=0) ;
        case 'clustering_func'
            v = v+1;
            clustering_func = input{v};
        case 'max2'
            v = v+1;
            max2 = input{v};
        case 'label'
            v = v+1;
            node_label = input{v};
        case 'nits'
            v = v+1;
            nits = input{v};
            otherwise
            error('pyramidalSGE:inputParser', 'Unsupported parameter: %s',input{v});
        end
        v = v+1;
    end
    if ~exist('max2', 'var')
        max2 = 7*ones(pyr_levels,1); % Max size of graphlets in terms of edges
    end
    max2 = uint32(max2);
end
