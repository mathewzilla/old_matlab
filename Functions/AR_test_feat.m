% AR_test_feat   Feature based classification
%
%   R = test_feature(data,coeffs)
%   Takes in a data file, and coefficients for an equation generated in 
%   train_feature.
%   The output is a radial distance estimate, rounded to the nearest R.

function class = AR_test_feat(data,coeffs)

  x = data;
  [peak ind] = max(x);
  
  % find width of peak
  down = 1;
  up = 0;
  m = 1;
  while x(ind+m)>0.05;
    down = down+1;
    m = m+1;
  end

  m = 1;
  while x(ind-m)>0.05;
    up = up +1;
    m = m+1;
  end
  
  width = down+up;
  coeffsR = coeffs(1,:);
  coeffsV = coeffs(2,:);

  radius = coeffsR(1)*(peak^4) + coeffsR(2)*(peak^3)*width + coeffsR(3)*(peak^3) + ...
           coeffsR(4)*(peak^2)*(width^2) + coeffsR(5)*(peak^2)*width + coeffsR(6)*(peak^2) ...
           + coeffsR(7)*peak*(width^3) + coeffsR(8)*peak*(width^2) + coeffsR(9)* ...
           peak*width + coeffsR(10)*peak + coeffsR(11)*(width^4) + coeffsR(12)*(width^3) ...
           + coeffsR(13)*(width^2) +coeffsR(14)*width + coeffsR(15);
  
  radius = round(radius);
  
  speed = coeffsV(1)*(peak^4) + coeffsV(2)*(peak^3)*width + coeffsV(3)*(peak^3) + ...
          coeffsV(4)*(peak^2)*(width^2) + coeffsV(5)*(peak^2)*width + coeffsV(6)*(peak^2) ...
          + coeffsV(7)*peak*(width^3) + coeffsV(8)*peak*(width^2) + coeffsV(9)* ...
          peak*width + coeffsV(10)*peak + coeffsV(11)*(width^4) + coeffsV(12)*(width^3) ...
          + coeffsV(13)*(width^2) +coeffsV(14)*width + coeffsV(15);
  
  speed = round(speed);
  
  % speed = ceil(speed);
  
  class = [radius;speed];