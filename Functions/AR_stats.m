% AR_stats takes the classification outputs of a classifier and computes
% the accuracy scores
function  [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class)

sz = size(class);
outR = zeros(1,length(indRadius));
outV = zeros(1,length(indRadius));
% Check if class output is actual R/V predictions
if max(max(indVelocity)) == max(max(class));
    outR = class(1,:);
    outV = class(2,:);
else
% if not use indRadius/Velocity to classify  
    for i = 1:length(indRadius);
        outR(i) = indRadius(1,class(1,i));
        if sz(1) == 2;   % if classifier makes a separate R/V prediction
            outV(i) = indVelocity(1,class(2,i));
        else
            outV(i) = indVelocity(1,class(1,i));
        end
    end
    
end

% SCALING 
if dset == 1;
    scR = 79; % Scaling by shortest R
    scRO = 1;  % Interval between Rs
    scS = 7.2; % Scaling by speed
    scSO = 7.2;
end
if dset == 2;
    scR = 1;
    scS = 1;
end
if dset == 3;
    scR = 1;
    scS = 1;
end

outR = outR*scO + scR;
indRadius = indRadius + scR;
outV = outV*scS;
indVelocity = indVelocity*scS;

errorR = outR - indRadius(1,:);
errorV = outV - indVelocity(1,:);

stat(1) = mean(abs(errorR)); % TO DO: Multiply by actual speed units
stat(2) = std(abs(errorR));
stat(3) = mean(abs(errorV));
stat(4) = std(abs(errorV));

