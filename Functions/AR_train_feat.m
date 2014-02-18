% AR_train_feat   Extract features and develop model from multivariate
% regression. Modified for publication in Evans et al 2014
% 'State of the art in whisker-based object localisation'.

%   Uses polyfitn.m which is included in the Functions folder, and was
%   written by John D'Errico as part of the sympoly toolbox
%
%   Coeffs = train_feature(DATA,index) for DATA an array.
%   Outputs coefficients of an equation that can be used to
%   classify new data.
%   index is the correct radial distance for each sample
%
%   Will include the degree of the regression (D) later. Default = 4 for now.

function out = AR_train_feat(data_train,indRadius,indVelocity,dset)

% First determine size of training set, then concatenate
% appropriately


indicesR = indRadius;
indicesV = indVelocity;
data_train = [data_train(3,:)];% ONLY USES ONE SET HERE. NEED TO REVISE
data_train = cell2mat(data_train);
[s,t] = size(data_train);

indicesR = indicesR(3,:);%,indicesR(2,:),indicesR(3,:)];
indicesV = indicesV(3,:);%,indicesV(2,:),indicesV(3,:)];
%% Extract features
features = zeros(t,2);

if dset == 1
    threshold = 0.05;
end

if dset == 2
    threshold = 0.05;
end

if dset == 3
    threshold = 4000;
end


% Find peak and its height
for i = 1: t;
    [peak ind] = max(data_train(:,i)); % used to be (:,i)
    features(i,1) = peak;
    %plot(data_train(:,i)); pause(0.3);
    
    % Find width of peak
    down = 1;
    up = 0;
    m = 1;
    while ind+m <= s && data_train(ind+m,i)>threshold;
        down = down+1;
        m = m+1;
    end
    
    m = 1;
    while ind-m >= 1 && data_train(ind-m,i)>threshold;
        up = up +1;
        m = m+1;
    end
    width = down+up;
    features(i,2) = width;
end

%% Run regression
% Bit of smoothing
if dset == 1;
    features(:,2) = smooth(features(:,2));
    
    % Sanity check the numbers - REVISE
    for i = 1:t; if features(i,2)<350; features(i,2) = features(i-1,2);end;end;
end

%% Generate polynomial linking peaks, widths, and the correct radial distance
p = polyfitn([features(:,1),features(:,2)],indicesR(:),4);
coeffsR = p.Coefficients;

%% Generate polynomial linking peaks, widths, and the correct speed
q = polyfitn([features(:,1),features(:,2)],indicesV(:),4);
coeffsV = q.Coefficients;

coeffs = [coeffsR;coeffsV];
out{1,1} = coeffs;
