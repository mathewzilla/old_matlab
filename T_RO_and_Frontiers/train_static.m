% train_static   Fit the static beam equation to the current
% experimental set up.
%
%   At the moment this is just a placeholder as I haven't worked
%   out a good way to do this. Right now it just generates the
%   equation coefficients



function eCoeffs = train_static(x,indRadius,indVelocity,theta);

% The way this will probably work is to take in appropriate
% parameters such as labelled training data, thetas of contact, and
% relevant whisker coefficients (wCoeffs), then fit the free
% parameters of the equation to the data.


% Whisker was 180mm long (ish), but it doesn't taper to
% nothing. The whisker is equivalent to a 229mm cone, truncated at
% 180mm.
% Contacts occur over 101mm range from the tip, outputs will
% therefore be from 80:180

L = 230;        %180; % 19/10/10 Need to work out effect of
                %truncated 229mm cone. Was 229
E = 1.5;        % 1.8 better than 2.3. WAS 2.4
rbase = 1.2;    % 0.9 get better results than 1. 1.1 'looks'
                % better, better sum of errors. Was 1
Ibase = pi.*(rbase.^4);
C = 3.*E.*(Ibase./4);
f = 0.002;     % Was 0.0033 (trying to fit to low V)

eCoeffs = [L E C f];
