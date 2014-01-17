% AR_run_temp script to run the template based classifer. Provides class
% predictions and clocked times

function [class,t_train,t_test] = AR_run_temp(data_train,data_test,method);
%% TRAINING
%   Old code
tic
fprintf(' ... Training ... ')
templates = AR_train_temp(data_train,method);
t_train = toc;



%% TESTING

% Old code
tic
fprintf('Testing ... \n')
% testing: template classifier  ----------------------------
clear class
for c = 1:length(data_test);
  class{c} = AR_test_temp(data_test{c},templates,method); 
  
end
class = cell2mat(class);
    

radius = class(1,:);
speed = class(2,:);

% --------------------------------------------------------
class(1,:) = radius;
class(2,:) = speed;

t_test = toc;