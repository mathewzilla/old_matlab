% train_feature   Extract features and develop model from multivariate regression.
%
%   Coeffs = train_feature(DATA,index) for DATA an array.
%   Outputs coefficients of an equation that can be used to
%   classify new data.
%   index is the correct radial distance for each sample
%
%   Will include the degree of the regression (D) later. Default = 4 for now.

function out = train_feature(x,indRadius,indVelocity,dims);

% First determine size of training set, then concatenate
% appropriately

nc1 = dims(1);
nc2 = dims(2);


indicesR = indRadius;
indicesV = indVelocity;
x = [x(2,:)];%,x(2,:),x(3,:)];
x = cell2mat(x);
[s,t] = size(x);

indicesR = [indicesR(2,:)];%,indicesR(2,:),indicesR(3,:)];
indicesV = [indicesV(2,:)];%,indicesV(2,:),indicesV(3,:)];
% Extract features
features = zeros(t,2);
% find peak and its height
for i = 1: t;
  [peak ind] = max(x(:,i));
  features(i,1) = peak;
  
  % find width of peak
  down = 1;
  up = 0;
  m = 1;
  while x(ind+m,i)>0.05;
    down = down+1;
    m = m+1;
  end

  m = 1;
  while x(ind-m,i)>0.05;
    up = up +1;
    m = m+1;
  end
  
  width = down+up;
  features(i,2) = width;
end

% Bit of smoothing
features(:,2) = smooth(features(:,2));

% Sanity check the numbers
for i = 1:t; if features(i,2)<350; features(i,2) = features(i-1,2);end;end;


 % generate polynomial linking peaks, widths, and the correct radial distance 
p = polyfitn([features(:,1),features(:,2)],indicesR(:),4);
coeffsR = p.Coefficients;

 % generate polynomial linking peaks, widths, and the correct speed
q = polyfitn([features(:,1),features(:,2)],indicesV(:),4);
coeffsV = q.Coefficients;

coeffs = [coeffsR;coeffsV];
out{1,1} = coeffs;

% Evaluating equation accross lots of data points
[X,Y] = meshgrid(linspace(min(features(:,1)),max(features(:,1)),100), ...
                 linspace(min(features(:,2)),max(features(:,2)),100));
        
modelR = polyvaln(p,[X(:),Y(:)]);

modelV = polyvaln(q,[X(:),Y(:)]);

%% UNCOMMENT TO SEE FEATURES (ONLY FOR COMPLETE SET) reshape and average
% USE polyn2sympoly(q) to generate polynomial

%features = reshape(features,nc1,nc2*3,2);
%feat1 = features(:,1:nc2,:);
%feat2 = features(:,nc2+1:2*nc2,:);
%feat3 = features(:,(2*nc2)+1:3*nc2,:);

%featAv = zeros(nc1,nc2,2);

%for i = 1:nc1; 
%  for j = 1:nc2; 
%    featAv(i,j,1) = mean([feat1(i,j,1),feat2(i,j,1),feat3(i,j,1)]);
%    featAv(i,j,2) = mean([feat1(i,j,2),feat2(i,j,2),feat3(i,j,2)]);
%  end;
%end;

%out{2,1} = featAv;



% Plotting;
%clf; hold all;
%plot3(features(:,1),features(:,2),indicesR(:),'.')
%plot3(X(:),Y(:),modelR(:),'+');



%  what the equation looks like (for classifying new data)
%  radius = coeffs(1)*(peak^4) + coeffs(2)*(peak^3)*width + coeffs(3)*(peak^3) + ...
%     coeffs(4)*(peak^2)*(width^2) + coeffs(5)*(peak^2)*width + coeffs(6)*(peak^2) ...
%     + coeffs(7)*peak*(width^3) + coeffs(8)*peak*(width^2) + coeffs(9)* ...
%     peak*width + coeffs(10)*peak + coeffs(11)*(width^4) + coeffs(12)*(width^3) ...
%     + coeffs(13)*(width^2) +coeffs(14)*width + coeffs(15);

% Next bit for playing with other data
%clf;
%for i = 1:4;
%  
% x = squeeze(featureBank(i,:,:,1)); % Peak deflection magnitude
% 
% x = reshape(x,1,2626); % reshape into vector
% 
% y = squeeze(featureBank(i,:,:,2)); % Deflection duration
% 
% y = reshape(y,1,2626); % reshape into vector
% 
% 
% gay = meshgrid(1:101,1:4:101); % generate array of radius labels
%
% gay = reshape(gay,1,2626); % reshape into vector
% 
% plot3(x,y,gay,'.')
% hold all;
% end
%
% p = polyfitn([x(:),y(:)],gay(:),4); % generate polynomial linking
%                                     % peaks, widths, and the
%                                     % correct radial distance 
%
% [X,Y] = meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100)) ;
% model = polyvaln(p,[X(:),Y(:)]);
%
%  
%plot3(X(:),Y(:),model(:),'+');