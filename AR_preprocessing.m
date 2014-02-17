function [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv)
% Preprocessing of datasets into appropriate training and test sets.

if dset == 1; % XY table data. 4 sets, 26 speeds, 101 radii
    
    % CODE FROM position_XY_feat.m in the first instance
    [ni,nc1,nc2] = size(data); [nt,nn] = size(data{1,1}); % 4,26,101 | 750,2
    
    ri = 1:ni;         % 1 to 4
    rc1 = 1:1:nc1;     % 1 to 26
    rc2 = 1:1:nc2;     % 1 to 101
    
    nc = nc1*nc2;      % total sample number per set - 2626
    rc = 1:nc;         % 1 to 2626
    subs = 1:3;        % subset for training set USE crossvalind INSTEAD
    
    %% onsets of contact. SEEMS REDUNDANT BUT KEEP IT AS IS USEFUL FOR OTHER SETS
    for c1 = rc1
        for c2 = rc2
            for i = ri
                t(i,c1,c2) = min(find(data{i,c1,c2}(:,1)>0.1))-100;
                
                
                
                %% OPTION: Take derivatives
                if deriv == 1;
                    data{i,c1,c2} = diff(data{i,c1,c2});
                    data{i,c1,c2}(nt,:) = data{i,c1,c2}(nt-1,:); % Fill in missing end value
                end
                
                if deriv == 2;
                    data{i,c1,c2} = diff(diff(data{i,c1,c2}));
                    data{i,c1,c2}(nt-1,:) = data{i,c1,c2}(nt-2,:); % Fill in missing end value
                    data{i,c1,c2}(nt,:) = data{i,c1,c2}(nt-1,:);
                end
                
                %% OPTION: Smooth
                if smoo == 1;
                    g = @(n,si) exp(-(-n:n).^2/(2*si^2))/sum( exp((-n:n).^2/(2*si^2)) );
                    data{i,c1,c2} = filter(g(50,15)/sum(g(50,15)),1,data{i,c1,c2});
                end
                
            end
        end
    end
    
    %% Separated files, not concatenated. ------
    
    % j determines which set is left out
    
    for j = subs;
        c= 0;
        for c1 = rc1
            for c2 = rc2
                c = c+1;
                rt = t(j+1,c1,c2)+(1:nt);
                data_train{j,c} = data{j+1,c1,c2}(rt,1);
                indRadius(j,c) = c2;
                indVelocity(j,c) = c1;
            end
        end
    end
    
    %% Concatenated files, if classifier needs it
    
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_train_c{c} = [];
            for k = subs
                rt = t(k+1,c1,c2)+(1:nt);
                data_train_c{c} = [data_train_c{c};data{k+1,c1,c2}(rt,1)]; % was vdata, derivative now happens in preprocessing
                
            end
        end
    end
    
    %% Test data: truncate
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_test{c} = [];
            for l = 1
                rt = t(l,c1,c2)+(1:nt);
                data_test{c} = [data_test{c};data{l,c1,c2}(rt,1)];
                test(c) = c;
            end
        end
    end
    
    
    
end



if dset == 2; %% Scratchbot (ROBIO expts). 3 speeds, 3 radii, 8 contacts
    % CODE FROM position_XY_feat.m in the first instance
    [nc1,nc2,ni] = size(data); [nt(1),nn] = size(data{1,1}); % 3,3,8 | 2001,2
    
    ri = 1:ni;         % 1 to 8 contacts
    rc1 = 1:1:nc1;     % 1 to 3 radii
    rc2 = 1:1:nc2;     % 1 to 3 speeds
    
    nc = nc1*nc2;      % total sample number per set - 9
    rc = 1:nc;         % 1 to 9
    subs = 1:7;        % subset for training set
    st = 1000;         % length of peak-aligned segment
    
    [nt(2),nn] = size(data{1,2});
    [nt(3),nn] = size(data{1,3});
    for c1 = rc1
        for c2 = rc2
            for i = ri
                t(c1,c2,i) = min(find(data{c1,c2,i}(:,1)>0.4))-100;
                % Pad with zeros as contacts are of different lengths
                if c2 >= 2;
                    data{c1,c2,i}(nt(c2):nt(1)-1,:)  = zeros(nt(1)-nt(c2),2);
                end
                
                %% OPTION: Take derivatives
                if deriv == 1;
                    data{c1,c2,i} = diff(data{c1,c2,i});
                    data{c1,c2,i}(nt(1)-1,:) = data{c1,c2,i}(nt(1)-2,:);% Fill in missing end value
                    data{c1,c2,i}(nt(1),:) = data{c1,c2,i}(nt(1)-1,:);
                end
                
                if deriv == 2;
                    data{c1,c2,i} = diff(diff(data{c1,c2,i}));
                    data{c1,c2,i}(nt(1)-2,:) = data{c1,c2,i}(nt(1)-3,:);% Fill in missing end value
                    data{c1,c2,i}(nt(1)-1,:) = data{c1,c2,i}(nt(1)-2,:);
                    data{c1,c2,i}(nt(1),:) = data{c1,c2,i}(nt(1)-1,:);
                end
                
                %% OPTION: Smooth
                if smoo == 1;
                    g = @(n,si) exp(-(-n:n).^2/(2*si^2))/sum( exp((-n:n).^2/(2*si^2)) );
                    data{c1,c2,i} = filter(g(50,15)/sum(g(50,15)),1,data{c1,c2,i});
                end
                
            end
        end
    end
    
    %% Separated files, not concatenated. ------
    
    % j determines which set is left out
    
    for j = subs;
        c= 0;
        for c1 = rc1
            for c2 = rc2
                c = c+1;
                rt = t(c1,c2,j+1)+(1:st);
                data_train{j,c} = data{c1,c2,j+1}(rt,1); % was (rt,1)
                indRadius(j,c) = c2;
                indVelocity(j,c) = c1;
            end
        end
    end
    
    %% Concatenated files, if classifier needs it
    
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_train_c{c} = [];
            for k = subs
                rt = t(c1,c2,k+1)+(1:st);
                data_train_c{c} = [data_train_c{c};data{c1,c2,k+1}(rt,1)]; % was (rt,1)
                
            end
        end
    end
    
    %% Test data: truncate
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_test{c} = [];
            for l = 1
                rt = t(c1,c2,l)+(1:st);
                data_test{c} = [data_test{c};data{c1,c2,l}(rt,1)]; % was (rt,1)
                test(c) = c;
            end
        end
    end
    
end



if dset == 3; %% Crunchbot data. 4 whiskers, 6 radii, 5 contacts
    % NEEDED TO COMBINE SET ACROSS WHISKERS. COULD LOOK AT INDIVIDUAL
    % WHISKERS LATER
    [nc1,nc2,ni] = size(data); [nt,nn] = size(data{1,1}); % 4,6,5 | 4000,2
    
    ri = 1:ni;         % 1 to 5 contacts
    rc1 = 1:1:nc1;     % 1 to 4 whiskers
    rc2 = 1:1:nc2;     % 1 to 6 radii
    
    nc = nc1*nc2*ni;   % total sample number per set - 120
    rc = 1:nc;         % 1 to 120
    subs = 1:7;        % subset for training set
    st = 1000;         % length of peak-aligned segment
    
    
    X  = data; clear data;
    %  for i = rc2
    data{i,:} = [X{1,i,:};
        
    
    for c1 = rc1
        for c2 = rc2
            for i = ri
                t(c1,c2,i) = min(find(data{c1,c2,i}(:,1)>0.4))-100;
                % Pad with zeros as contacts are of different lengths
                if c2 >= 2;
                    data{c1,c2,i}(nt(c2):nt(1)-1,:)  = zeros(nt(1)-nt(c2),2);
                end

                %% OPTION: Take derivatives
                if deriv == 1;
                    data{c1,c2,i} = diff(data{c1,c2,i});
                    data{c1,c2,i}(nt(1)-1,:) = data{c1,c2,i}(nt(1)-2,:);% Fill in missing end value
                    data{c1,c2,i}(nt(1),:) = data{c1,c2,i}(nt(1)-1,:);
                end

                if deriv == 2;
                    data{c1,c2,i} = diff(diff(data{c1,c2,i}));
                    data{c1,c2,i}(nt(1)-2,:) = data{c1,c2,i}(nt(1)-3,:);% Fill in missing end value
                    data{c1,c2,i}(nt(1)-1,:) = data{c1,c2,i}(nt(1)-2,:);
                    data{c1,c2,i}(nt(1),:) = data{c1,c2,i}(nt(1)-1,:);
                end

                %% OPTION: Smooth
                if smoo == 1;
                    g = @(n,si) exp(-(-n:n).^2/(2*si^2))/sum( exp((-n:n).^2/(2*si^2)) );
                    data{c1,c2,i} = filter(g(50,15)/sum(g(50,15)),1,data{c1,c2,i});
                end

            end
        end
    end
    
    %% Separated files, not concatenated. ------
    
    % j determines which set is left out
    
    for j = subs;
        c= 0;
        for c1 = rc1
            for c2 = rc2
                c = c+1;
                rt = t(c1,c2,j+1)+(1:st);
                data_train{j,c} = data{c1,c2,j+1}(rt,1); % was (rt,1)
                indRadius(j,c) = c2;
                indVelocity(j,c) = c1;
            end
        end
    end
    
    %% Concatenated files, if classifier needs it
    
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_train_c{c} = [];
            for k = subs
                rt = t(c1,c2,k+1)+(1:st);
                data_train_c{c} = [data_train_c{c};data{c1,c2,k+1}(rt,1)]; % was (rt,1)
                
            end
        end
    end
    
    %% Test data: truncate
    c = 0;
    for c1 = rc1
        for c2 = rc2
            c = c + 1; data_test{c} = [];
            for l = 1
                rt = t(c1,c2,l)+(1:st);
                data_test{c} = [data_test{c};data{c1,c2,l}(rt,1)]; % was (rt,1)
                test(c) = c;
            end
        end
    end
    
    
    
    


    
end



% if dset == 4; %% Extra Scratchbot dataset
%
%
%
%
%
%
%
% end


%return