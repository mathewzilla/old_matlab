% AR_run_static script to run the static beam equation based classifier. Provides class
% predictions and clocked times

function [class,t_train,t_test] = AR_run_static(data_train,data_test,indRadius,indVelocity);

%% TRAINING
tic
fprintf(' ... Training ... ')
%% Generate protraction angle (theta) for each R-S pair
theta = AR_ProtractionAngle;

%% Set parameter values for static beam equation
L = 180; rbase = 1; E = 2.4; % 180, 1. 2.4/ 230, 1.2, 1.5
Ibase = pi.*(rbase.^4);
C = 3.*E.*(Ibase./4);
f = 0.0031; %0.0001:0.0001:0.01; %0.0005: 0.0005: 0.005;     % Was
%0.0031 / 0.002

eCoeffs = [L E C f];
t_train = toc;

%% TESTING

% Old code
tic
fprintf('Testing ... \n')
% testing: static beam equation based classifier  ----------------------------
clear class
for c = 1:length(data_test);
    class{c} = AR_test_static(data_test{c},theta(c),eCoeffs);
    
end
class = cell2mat(class);
% Normalising the output
radius = max(class)-class;    % Output is wrong way round initially
radius = radius./max(radius); % Normalise to 1
radius = (radius*100) +1;     % Scale up to 1-101 range
radius = round(radius);       % Round to nearest integer (for comparison with other methods)

class(1,:) = radius;
class(2,:) = 13*ones(1,2626); % NO SPEED OUTPUT FOR THIS CLASSIFIER

t_test = toc;



%% LOOP TO FIND BEST PARAMETERS - SHOULD OPTIMIZE PROPERLY. IMPORTANT FOR COMPUTING CLASSIFIER TIMING
% pish = zeros(4,15*26*11*10);
% meh =  zeros(3,15*26*11*10);
%  for L = 170 :5: 240;        %180; % 19/10/10 Need to work out effect of
%                %truncated 229mm cone. Was 229
%  for E = 0.5 :0.1:3 ;        % 1.8 better than 2.3. WAS 2.4
%    for rbase = 0.5:0.1:1.5;    % 0.9 get better results than 1. 1.1 'looks'
%                % better, better sum of errors. Was 1
%      Ibase = pi.*(rbase.^4);
%      C = 3.*E.*(Ibase./4);
% THEN TEST THESE VALUES (not done here)
% end; end; end;