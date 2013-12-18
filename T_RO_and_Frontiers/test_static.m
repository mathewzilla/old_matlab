% test_static.m Static beam equation based classification
%
%   R = test_static(data,theta)
%   Takes in a data file, and whisk angle (theta). 
%   On the XY table this has been calculated from contact speeds, 
%   locations and times. On a whisking robot it will be
%   controlled. 
%   Otherwise it must be guessed.
%   The output is a radial distance estimate, rounded to the
%   nearest R.

function class = test_static(data,theta,eCoeffs);

  x = data;
  % [peak ind] = max(x);

  
  L = eCoeffs(1); %229;        %180; % 19/10/10 Need to work out effect of truncated 229mm cone
  E = eCoeffs(2); %2.3;        % 1.8 better than 2.3
  %rbase = params(3); %1.0;    % 0.9 get better results than 1. 1.1 'looks' better, better sum of errors
  %Ibase = pi.*(rbase.^4);
  C = eCoeffs(3);
  f = eCoeffs(4); %0.0033;

  F = f*(downsample(squeeze(x),10)); %0.0033 is the best

  M = max(F);
  % NEED TO REMEMBER THAT THETA MUST BE CALCULATED TOO..
  %theta = Theta(S,R);

  %thetaDot = diff(theta);
  %Mdot = diff(M);
  %radius(R) = R + 79;
  D = (C.*theta.*L)./((C.*theta)+(M.*L));
  
  class = D;
  