%   AR_test_temp   Template based classification

%   R = test_feature(data,coeffs)
%   Takes in a data file, and trained templates. Probably want to add an
%   option for method
%   The output is a radial distance estimate, rounded to the nearest R.

function class = AR_test_temp(data_test,templates,method); 

% Included option for template 'method'
% 1 = low pass filtered
% 2 = average template
% 3 = Velocity
% 4 = Velocity of single trial
% 5 = Frequency template

%% format template of test data, then find sum squared error
x = data_test; nx = length(x); rx = 1:nx; 
[nc ng] = size(templates); rc = 1:nc;  

%% Low pass filter
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

%% Average template
if method == 2;      
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

%% First derivative of average template REDUNDANT?
if method == 3;      
  temp = diff(x);

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in;
end

%% First derivative of single trial REDUNDANT?
if method == 4;      
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

%% Frequency spectrum template
if method == 5;
  temp =  abs(real(fft(x(501:750))));
  temp = smooth(smooth(temp));
  temp(1) = 0;
  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  class = in; %indRadius(in);
end