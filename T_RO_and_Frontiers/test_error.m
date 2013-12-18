% test_error   error based classification
%
%   error = test_error(data,features)
%   Takes in a data file, and features from train_feature
%   The output is an array of errors, representing the difference
%   between the input features and those in the array.

function error = test_error(data,features);

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

for i = 1:26;
  for j = 1:26;
    error(i,j,1) = abs(peak-features(i,j,1));
    error(i,j,2) = abs(width-features(i,j,2));
  end;
end;

error = reshape(error,676,2);