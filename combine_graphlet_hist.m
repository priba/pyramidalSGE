function [ hist ] = combine_graphlet_hist( histogram, t, combine_graphlet )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if strcmpi(combine_graphlet, 'combine')
        hist = [histogram{1:t}] ;
    else
        hist = histogram{t} ;
    end ;

end

