% train_template	Store templates for later
% classification. Just stock as of 13/5/11, but will want to do
% something extra.
%
%   out = train_template(DATA) for DATA an array.
%   Outputs an array of templates. Might do something clever with
%   smoothing, phase diagrams or even 'dynamic time warping'!

% Adaptive 'EM' version maybe a bit wooly.
%

function templates = train_template(data,method,smoo); % No other parameters
                                          % needed in this
                                          % version. Maybe want
                                          % smoothing etc later

                                          % Include option for
                                          % template 'method' 
                                          % 1 = low pass filtered
                                          % 2 = average template
                                          % 3 = DTW
                                          % 4 = Velocity
                                          % 5 = Velocity of single
                                          % trial 
                                          % 6 = Corr Coeff
                                          % 7 = Frequency template

[nt nn] = size(data); %[nd nnd] = size(d);


%if nnd==1; d = repmat(d,[1,nn]); end

if method == 1;
  f = 1000;                % Sample rate. Should probably have this
                           % as input

  % Low pass filter the samples
    fNorm = 20/(f/2);        % Low pass cutoff at 20Hz. Should input
                           % this too
  [b a] = butter(10, fNorm, 'low');

  for i = 1:nn;
 %   clear x;
    clear y;
    for j = 2;  %:nt;
   %   y(1,:) = data{j,i};
      y(1,:) = filter(b,a, data{j,i});
    end
    if smoo;
      y = smooth(y);
    end
    
    templates(i,:) = y; 
    
  end
end

if method == 2;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 1:3;  %:nt;
      x(j,:) = data{j,i};
    end
 
    for k = 1:length(x);
     y(k) = mean(x(:,k));
    end
        
    if smoo;
      y = smooth(y);
    end
    templates(i,:) = y; 
    
  end
end

if method == 3;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 1;  %:nt;

    end
    %for k = 1:length(y);
    % z(k) = mean(y(:,k));
    %end
    templates(i,:) = y; % was z
    
  end
end

if method == 4;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 1:3;  %:nt;
      x(j,:) = data{j,i};
    end
 
    for k = 1:length(x);
     y(k) = mean(x(:,k));
    end
    
    y = diff(y);
        
    if smoo;
      y = smooth(y);
    end
    templates(i,:) = y; 
    
  end
end

if method == 5;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 2;  %:nt;
      x = data{j,i};
    end
    
    y = diff(x);
        
    if smoo;
      y = smooth(y);
    end
    templates(i,:) = y; 
    
  end
end

if method == 6;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 2;  %:nt;
      x = data{j,i};
    end
    
    if smoo;
      x = smooth(x);
    end
    templates(i,:) = x; 
    
  end
end

if method == 7;

  for i = 1:nn;
    clear x;
    clear y;
    for j = 2;
      x = abs(real(fft(data{j,i}(501:750))));
    end
    
      x(1) = 0;
      if smoo;
        x = smooth(smooth(x));
      end
      templates(i,:) = x; 
  end
end
  
