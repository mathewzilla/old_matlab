% AR_run_SNB script to run the Naive Bayes classifer. Provides class
% predictions and clocked times

function [class,t_train,t_test] = AR_run_SNB(data_train_c,data_test);
%% TRAINING

%   TO DO:  AR_train('SNB'); New universal code?
%   Old code
d = linspace(-0.1,0.1,501); d = d(:); % set params for histogram, then rotate
ns = 10;                              % smooth over 10 samples
clear logl;
tic
fprintf(' ... Training ... ')
for c = 1:length(data_train_c);
    logl{c} = AR_train_SNB(data_train_c{c},d,ns);
end
t_train = toc;

%% TESTING

% Old code
clear logp;
tic
fprintf('Testing ... \n')
for c = 1:length(data_train_c);
    for c1 = 1:length(data_train_c);
        logp{c}(c1) = AR_test_SNB(data_test{c},d,logl{c1});
    end
    [ig,class(c)] = max(squeeze(logp{c}),[],2);
end
t_test = toc;