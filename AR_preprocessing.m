function [data_train,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset)
% Preprocessing of datasets into appropriate training and test sets

if dset == 1; % XY table data. 4 sets, 26 speeds, 101 radii
    
    % CODE FROM position_XY_feat.m in the first instance
    [ni,nc1,nc2] = size(data); [nt,nn] = size(data{1,1}); % 4 by 26 by
    % 101. 750 by 2
    ri = 1:ni; rc1 = 1:1:nc1; rc2 = 1:1:nc2; % was every 4th trial (radius)
    nc1 = length(rc1); nc2 = length(rc2); nc = nc1*nc2; rc = 1:nc;
    
    % onsets
    for c1 = rc1
        for c2 = rc2
            for i = ri
                t(i,c1,c2) = min(find(data{i,c1,c2}(:,1)>0.1))-100;
            end
        end
    end
    
    
    
    % Separated files, not concatenated. ------
    
    for i = 1:3;
        c= 0;
        for c1 = rc1
            for c2 = rc2
                c = c+1;
                rt = t(i+1,c1,c2)+(1:750);
                data_train{i,c} = data{i+1,c1,c2}(rt,1);
                indRadius(i,c) = c2;
                indVelocity(i,c) = c1;
            end
        end
    end
    
    % test data: truncate
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_test{c} = [];
            for i = 1
                rt = t(i,c1,c2)+(1:750);
                data_test{c} = [data_test{c};data{i,c1,c2}(rt,1)];
                test(c) = c;
            end
        end
    end
    
    
    
end



if dset == 2; % Scratchbot (ROBIO expts). 3 speeds, 3 radii, 8 contacts
    
    
    
    
    
    
    
end



if dset == 3; % Crunchbot data. 4 whiskers, 6 radii, 5 contacts
    
    
    
    
    
    
    
end



if dset == 4; % XY table data
    
    
    
    
    
    
    
end


return