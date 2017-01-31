function [  ] = classify_dataset( dataset_name, varargin )
%CLASSIFY_DATASET Summary of this function goes here
%   Detailed explanation goes here

    %% Parameters
        
    % Check hostname to know the file structure for adding the packages
    [~, name] = system('hostname');
    
    if ismember(lower(name), {'cvc206'} ) 
        % Dataset path
        params.p_data = '/home/adutta/Workspace/Datasets/STDGraphs' ;

        % External libreries (considering them to be on the userpath folder
        params.libraries = {
            '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools/matlab_bgl',...
            '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools/libsvm/matlab',...
            '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools/vlfeat'
            } ;

             % Project folders
        params.folders = { 'clustering',...
            '/home/adutta/Dropbox/Personal/Workspace/AdditionalTools/random_graphlet1'
            } ;
    else
        % Dataset path
        params.p_data = [ pwd filesep 'dataset' filesep ] ;

        % External libreries (considering them to be on the userpath folder
        params.libraries = { 'matlab_bgl', 'libsvm/matlab', 'vlfeat' } ;

        % Project folders
        params.folders = { 'clustering', 'random_graphlet1' } ;
    end ;
    
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
    % Graph sizes (in terms of edges)
    params.MAX2 = uint32(8);

    %% Addpaths
    folders_paths = cellfun(@(x) strcat(pwd, filesep,x),params.folders, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath([folders_paths{:} ])
    
    user_path = strrep(userpath,';',''); % Windows
    user_path = strrep(user_path,':',''); % Linux
    
    folders_paths = cellfun(@(x) strcat(user_path, filesep,x),params.libraries, 'UniformOutput', false ) ;
    folders_paths = cellfun(@genpath,folders_paths, 'UniformOutput', false ) ;
    addpath( [folders_paths{:} ])

    clear user_path folders_paths ;
    
    %% Default values
    [epsi, del, pyr_levels, pyr_reduction, edge_thresh, clustering_func, nits, VERBOSE]...
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

    M = uint32(ceil(2*(params.T(1:params.MAX2)*log(2)+log(1/del))/epsi^2));
    
    % Initialize storage
    global_var.hash_codes_uniq = cell(params.MAX2-2,1);
    global_var.idx_image = cell(params.MAX2-2,1);
    global_var.idx_bin = cell(params.MAX2-2,1);
    
    %% Iterate whole dataset
    for i = 1:ngraphs
        if VERBOSE
            fprintf('Graph: %d.\n',i);
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
            [ global_var ] = graphlet_embedding(pyr_graph{j} , i , M , global_var, params ) ;
        end ;
        
        
        if VERBOSE
            toc
        end ;
    end;
    
    dim_hists = cellfun(@(x) size(x,1) ,global_var.hash_codes_uniq);
    clear global_var.hash_codes_uniq;
    
    %% Compute histograms and kernels
    histograms = cell(1,params.MAX2-2);

    KM_train = zeros(ngraphs,ngraphs,params.MAX2-2);
    KM_test = zeros(ngraphs,ngraphs,params.MAX2-2);

    for i = 1:params.MAX2-2

        histograms{i} = sparse(global_var.idx_image{i},global_var.idx_bin{i},1,ngraphs,dim_hists(i));
        histograms{i} = bsxfun(@times, histograms{i},1./(sum(histograms{i},2)+eps));

        X_train = histograms{i};
        X_test = histograms{i};

        KM_train(:,:,i) = vl_alldist2(X_train',X_train','KL1');
        KM_test(:,:,i) = vl_alldist2(X_test',X_train','KL1');

    end;
    
    clear global_var;
    
    %% Evaluate
    % Evaluate nits times to get the accuracy mean and standard deviation
    accs = zeros(nits,params.MAX2-2);
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

        for i = 1:params.MAX2-2

            w = ones(1,i);
            w = w/sum(w);

            K_train = [(1:ntrain_set)' KM_train(train_set,train_set,i)];
            K_test = [(1:ntest_set)' KM_test(test_set,train_set,i)];

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
            accs(it,i) = acc(1);

        end;
    end ;
    
    % Mean and standard deviation
    maccs = mean(accs);
    mstds = std(accs)./sqrt(nits);
    
    % Save results
    fileID = fopen([params.out dataset_name , '.txt'],'a') ;
    fprintf(fileID,params.headerSpec, dataset_name, ngraphs, nclasses, nits) ;
    fprintf(fileID,params.sgeSpec, epsi , del) ;
    fprintf(fileID,params.pyrSpec, pyr_levels , pyr_reduction , edge_thresh , func2str(clustering_func)) ;
    for i = 1:params.MAX2-2
        fprintf(fileID, 't = %d \t %.2f \\pm %.2f \n',i+2, maccs(i),mstds(i));
        if VERBOSE
            fprintf('t = %d \t %.2f\\pm %.2f \n', i+2 , maccs(i),mstds(i));
        end
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
        otherwise
            error('pyramidalSGE:incorrectDataset',...
                'Error. \nDatabase %s not accepted.', dataset_name)
    end;
end

%% Default values
function [eps, del, pyr_levels, pyr_reduction, edge_tresh, clustering_func,...
    nits, VERBOSE] = input_parser( input )
    VERBOSE = 0 ;
    eps = 0.1 ;
    del = 0.1 ;
    pyr_levels = 1 ;
    pyr_reduction = 2 ;
    edge_tresh = 0 ;
    clustering_func = @girvan_newman ;
    nits = 10 ;
    
    % Parse optional input parameters
    v = 1;
    while v < numel(input)
        switch input{v}
        case 'VERBOSE'
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
        case 'nits'
            v = v+1;
            nits = input{v};
            otherwise
            error('pyramidalSGE:inputParser', 'Unsupported parameter: %s',input{v});
        end
        v = v+1;
    end
end
