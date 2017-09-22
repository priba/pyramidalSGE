%% Original graph
% Random seed
rng(0);

% Create a random graph
num_nodes = 1000;
prob_edge = 0.001;
A = rand(num_nodes)<prob_edge;
[I,J] = find(A);

% List of edges
L = [uint32(I),uint32(J)]; % Adjacency edges
L_h = []; % Hierarchical edges

% Parameters
num_graphlets = uint32(10000);
max_len_graphlets = uint32(3);

% Generate random graphs
tic;
graphlets = generate_random_hierarchical_graphlets(L,L_h,num_graphlets,max_len_graphlets);
toc;

%% Hierarchical graph
H1 = A; % Original graph
H12 = rand(num_nodes, round(num_nodes/2))<prob_edge; % Hierarchical edges
H2 = rand(round(num_nodes/2))<prob_edge; % 1st Abstract level
H = [H1, H12; H12', H2]; % Whole adjacency matrix of the hierarchical representation

[I1,J1] = find(H1);
[I12,J12] = find(H12); J12 = J12 + size(H1,1);
[I2,J2] = find(H2);

% List of edges
L = [uint32(I1),uint32(J1); uint32(I2),uint32(J2)];  % Adjacency edges
L_h = [uint32(I12),uint32(J12)]; % Hierarchical edges

% Generate random graphs
tic;
graphlets = generate_random_hierarchical_graphlets(L,L_h,num_graphlets,max_len_graphlets);
toc;
