%   AR_test_temp   Template based classification

%   R = test_feature(data,coeffs)
%   Takes in a data file, and trained templates. Probably want to add an
%   option for method
%   The output is a radial distance estimate, rounded to the nearest R.

function class = AR_test_temp(data_test,templates,method); 

%% format template of test data, then find sum squared error
x = data_test; nx = length(x); rx = 1:nx; 
[nc ng] = size(templates); rc = 1:nc;  

% format templates of test
% Filter as for training set
if method == 1;            % LOW PASS FILTER
  f = 1000;                % Sample rate. 
  fNorm = 20/(f/2);        % Low pass cutoff at 20Hz.
  [b a] = butter(10, fNorm, 'low');
  temp = filter(b,a,x);
%  temp = x;

  
  % template classifier using RMS error
  [nt nn] = size(x);
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in; 
  % class = indRadius(in);
end

if method == 2;      % AVERAGE TEMPLATE
  temp = x;

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in; 
  % class = indRadius(in);
end

if method == 3;      % DTW (DIDN"T DO IT IN THE END)
  temp = x;

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in; 
  % class = indRadius(in);
end

if method == 4;      % VELOCITY FROM AVERAGE TRIAL
  temp = diff(x);

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = indRadius(in);
end

if method == 5;      % VELOCITY OF A SINGLE TRIAL
  temp = diff(x);
  
  % template classifier using RMS error
  [nt nn] = size(x);
  for j = rc;         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in; 
  % class = indRadius(in);
end

if method == 6;     % CORRELTAION COEFFICIENT
  temp = x;

  % template classifier using corr coeff
  [nt nn] = size(x); 
  for j = rc;       % for 1:676 (26*26)
    Cor(j) = corr(temp,templates(j,:)'); 
  end; 
  [k in] = max(Cor); 
  class = in; 
  % class = indRadius(in); 

end

if method == 7;
  temp =  abs(real(fft(x(501:750))));
  
  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = in; %indRadius(in);
end