% AR_ProtractionAngle
% Script to work out the angle and protraction distance (opposite side of
% the triangle to the angle)          \--D--|
%                                      \    |
%                                     H \   | R
%                                        \~~|
%                                         \A|
%                                          \|
%                                         =====

function [theta,eCoeffs] = AR_ProtractionAngle(dset)
if dset == 1;
    % Set up parameter ranges for feeding through the trig.
    R = linspace(80,180,101); % All mm along the whisker
    S = linspace(36,216,26);  % All speeds in the range
    D = S * 0.3;              % Distance of the robot, given protraction lasts 0.3 seconds
    H = zeros(26,101);        % Hypotenuse length for each S-R pair
    A = zeros(26,101);        % Angle for each S-R pair
    Theta = zeros(26,101);
    theta = zeros(1,2626);    % One long vector to match data_test
    k = 0;
    for i = 1:26;
        for j = 1:101;
            k = k+1;
            H(i,j) = sqrt(D(i).^2 + R(j).^2);
            A(i,j) = asin(D(i)./H(i,j));
            Theta(i,j) = D(i)./R(j);
            theta(k) = Theta(i,j);
        end
    end
    %% Set parameter values for static beam equation
    L = 185; rbase = 1; E = 2.4; % ABS 180, 1. 2.4/ 230, 1.2, 1.5
    Ibase = pi.*(rbase.^4);
    C = 3.*E.*(Ibase./4);
    f = 0.0031; %0.0001:0.0001:0.01; %0.0005: 0.0005: 0.005;     % Was
    %0.0031 / 0.002
    
    eCoeffs = [L E C f];
end

if dset == 2;
    % Set up parameter ranges for feeding through the trig.
    R = [70 100 130]; % in mm along the whisker
    S = [2 4 6];     % Whisk speed in Hz
    theta = deg2rad(15)*ones(1,9); % From Evans et al Frontiers 2013
    %% Set parameter values for static beam equation
    L = 185; rbase = 1; E = 25; % Fiberglass
    Ibase = pi.*(rbase.^4);
    C = 3.*E.*(Ibase./4);
    f = 0.1;    % Back of the envelope stuff. Need to callibrate properly with data
    
    eCoeffs = [L E C f];
end

if dset == 3;
    % Set up parameter ranges for feeding through the trig.
    R = [90 100 110 120 130 140]; % in mm along the whisker
    S = 1:4;     % whisker identity (not used)
    theta = deg2rad(15)*ones(1,24); % From Evans et al Frontiers 2013
    %% Set parameter values for static beam equation
    L = 160; rbase = 1.45; E = 4.89; % Fiberglass
    Ibase = pi.*(rbase.^4);
    C = 3.*E.*(Ibase./4);
    f = 0.0001;    % Back of the envelope stuff. Need to callibrate properly with data
    
    eCoeffs = [L E C f];
end