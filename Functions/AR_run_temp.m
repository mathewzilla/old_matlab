% AR_run_temp script to run the template based classifer. Provides class
% predictions and clocked times

% Included option for template 'method'
% 1 = low pass filtered
% 2 = average template
% 3 = Velocity
% 4 = Velocity of single trial
% 5 = Frequency template


function [class,t_train,t_test] = AR_run_temp(data_train,data_test,indRadius,indVelocity,method);


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
    

radius = indRadius(1,class(:));
speed = indVelocity(1,class(:));

% --------------------------------------------------------
class(1,:) = radius;
class(2,:) = speed;

t_test = toc;