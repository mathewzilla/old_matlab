% Script to work out the angle and protraction distance (opposite side of
% the triangle to the angle)          \--D--|
%                                      \    |
%                                     H \   | R
%                                        \~~|
%                                         \A|
%                                          \|
%                                         =====


R = linspace(80,180,101);
S = linspace(36,216,26);
D = S * 0.3;
H = zeros(26,101);
A = zeros(26,101); 
AngleRad = zeros(26,101); 
AngleDeg = zeros(26,101);
Theta = zeros(26,101);
for i = 1:26;
    for j = 1:101;
        H(i,j) = sqrt(D(i).^2 + R(j).^2);
        A(i,j) = asin(D(i)./H(i,j));
        Theta(i,j) = D(i)./R(j);
    end
end

yD = zeros(26,300);
ThetaVec = zeros(26,101,300);
% Generate vector version for 
for i = 1:26;
    for j = 1:101;
        for k = 1:300;
            yD(i,k) = (D(i)./300).*k;
            ThetaVec(i,j,k) = yD(i,k)./R(j);
        end
    end
end


AngleRad  = A;
AngleDeg = AngleRad.*180./pi;

% To see the point where the hall sensor should be maxed out
for i = 1:26; 
    for j = 1:101; 
        if AngleDeg(i,j) >= 17 && AngleDeg(i,j) <=18; 
            AngleDeg(i,j) = 40; 
        end; 
    end;
end

% Load the data for plotting
%load /Users/matevans/Documents/Matlab_code/SAB2010/dataSetAll.mat
%load /Users/matevans/Documents/Matlab_code/SAB2010/featuresHomeDown/Normfeats

% figure(5);
% subplot(2,1,1); 
% imagesc(AngleDeg), axis image, set(gca,'TickDir','out');
% subplot(2,1,2);
% imagesc(normFeats(:,:,1)), axis image, set(gca,'Xdir','reverse','TickDir','out');


