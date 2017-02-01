function [ combs ] = allcomb( v1 , v2 )
%ALLCOMB Summary of this function goes here
%   http://stackoverflow.com/questions/21895335/generate-a-matrix-containing-all-combinations-of-elements-taken-from-n-vectors

    vectors = { v1, v2 }; %// input data: cell array of vectors

    n = numel(vectors); %// number of vectors
    combs = cell(1,n); %// pre-define to generate comma-separated list
    [combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); %// the reverse order in these two
    %// comma-separated lists is needed to produce the rows of the result matrix in
    %// lexicographical order 
    combs = cat(n+1, combs{:}); %// concat the n n-dim arrays along dimension n+1
    combs = reshape(combs,[],n); %// reshape to obtain desired matrix
end
