% AR_run_feat script to run the feature based classifer. Provides class
% predictions and clocked times

function [class,t_train,t_test] = AR_run_feat(data_train,data_test,indRadius,indVelocity,dset);
%% TRAINING

%   Old code

tic
fprintf(' ... Training ... ')
out = AR_train_feat(data_train,indRadius,indVelocity,dset);
t_train = toc;

coeffs = out{1,1};


%% TESTING

% Old code
tic
fprintf('Testing ... \n')
% testing: feature classifier  ----------------------------
clear class
for c = 1:length(data_test);
    class{c} = AR_test_feat(data_test{c},coeffs);
end
class = cell2mat(class);

radius = class(1,:);
speed = class(2,:);

%% Bound max and min of reports to reasonable values
for i = 1:length(data_test);
    if radius(i)<1;
        radius(i) = 1;
    end;
    if radius(i)>101;
        radius(i) = 101;
    end;
end;

for i = 1:length(data_test);
    if speed(i)<1;
        speed(i) = 1;
    end;
    if speed(i)>26;
        speed(i) = 26;
    end;
end;
% --------------------------------------------------------
class(1,:) = radius;
class(2,:) = speed; % ONLY ONE THAT MAKES AN EXPLICIT PREDICTION OF SPEED

t_test = toc;