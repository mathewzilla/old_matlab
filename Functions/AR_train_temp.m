% AR_train_temp    Store templates for later classification. Modified for 
%   publication in Evans et al 2014
%   'State of the art in whisker-based object localisation'.

%   out = AR_train_temp(data_train) for DATA an array.
%   Outputs an array of templates. Not very clever! This is basically for 
%   consistency and timing.
%   Old code executed all kinds of preprocessing, but now that is handled 
%   by AR_preprocessing.m
%

function out = AR_train_feat(data_train)


out = data_train;