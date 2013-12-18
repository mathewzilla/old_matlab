% test_template   template based classification
%
%   R = test_feature(data,coeffs)
%   Takes in a data file, and trained templates. Probably want to add an
%   option for method
%   The output is a radial distance estimate, rounded to the nearest R.


% JUST CODE FROM CLASSIFY_TEMPLATE.M AT THE MOMENT 13/5/11
function c = test_template(testdata,templates,indRadius,method,smoo); % Add option for template
                                        % method later. e.g EM,
                                        % DTW, Phase plane etc


%% format template of test data, then find sum squared error
x = testdata; nx = length(x); rx = 1:nx; 
[nc ng] = size(templates); rc = 1:nc;  

% format templates of test
% Filter as for training set
if method == 1;            % LOW PASS FILTER
  f = 1000;                % Sample rate. 
  fNorm = 20/(f/2);        % Low pass cutoff at 20Hz.
  [b a] = butter(10, fNorm, 'low');
  temp = filter(b,a,x);
%  temp = x;
  if smoo;
    temp = smooth(temp);
  end
  
  % template classifier using RMS error
  [nt nn] = size(x);
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = in; 
  % c = indRadius(in);
end

if method == 2;      % AVERAGE TEMPLATE
  temp = x;
  
  if smoo;
    temp = smooth(temp);
  end

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = in; 
  % c = indRadius(in);
end

if method == 3;      % DTW (DIDN"T DO IT IN THE END)
  temp = x;

  if smoo;
    temp = smooth(temp);
  end

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = in; 
  % c = indRadius(in);
end

if method == 4;      % VELOCITY FROM AVERAGE TRIAL
  temp = diff(x);

  if smoo;
    temp = smooth(temp);
  end

  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = indRadius(in);
end

if method == 5;      % VELOCITY OF A SINGLE TRIAL
  temp = diff(x);

  if smoo;
    temp = smooth(temp);
  end

  % template classifier using RMS error
  [nt nn] = size(x);
  for j = rc;         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = in; 
  % c = indRadius(in);
end

if method == 6;     % CORRELTAION COEFFICIENT
  temp = x;

  if smoo;
    temp = smooth(temp);
  end

  % template classifier using corr coeff
  [nt nn] = size(x); 
  for j = rc;       % for 1:676 (26*26)
    Cor(j) = corr(temp,templates(j,:)'); 
  end; 
  [k in] = max(Cor); 
  c = in; 
  % c = indRadius(in); 

end

if method == 7;
  temp =  abs(real(fft(x(501:750))));
  
  if smoo;
    temp = smooth(smooth(temp));
  end
  
  % template classifier using RMS error
  [nt nn] = size(x); 
  for j = rc         % for 1:676 (26*26)
    Rsq(j) = sum(abs(templates(j,:) - temp'));
  end

  [ig,in] = min(Rsq); 
  c = indRadius(in);
end