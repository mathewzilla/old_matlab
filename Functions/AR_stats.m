% AR_stats takes the classification outputs of a classifier and computes
% the accuracy scores
function  [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class)

sz = size(class);
outR = zeros(1,length(indRadius));
outV = zeros(1,length(indRadius));
for i = 1:length(indRadius);
    outR(i) = indRadius(1,class(1,i));
    if sz(1) == 2;
    outV(i) = indVelocity(1,class(2,i));
    else
        outV(i) = indVelocity(1,class(1,i));
    end
end
errorR = outR - indRadius(1,:);
errorV = outV - indVelocity(1,:);

stat(1) = mean(abs(errorV)); % TO DO: Multiply by actual speed units
stat(2) = std(abs(errorV));
stat(3) = mean(abs(errorR));
stat(4) = std(abs(errorR));

