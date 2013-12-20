% This code runs the old classifiers on the 3 different robot datasets,
% under all the relevant conditions for publication in Evans et al 2014
% 'State of the art in whisker-based object localisation'.
% Puts all of the outputs in an array for easy copying into the paper.
% AR_figs.m generates the publcation quality figures. See README.txt for
% more details.

clear all; clc;

% load the three (not four) data sets, one at a time.
for dset = 1;% :3; other numbers don't work yet
    clear data
    if dset == 1;load('Data/data_XY','data');   % XY table
    end;                                        % 4 sets, 26 speeds, 101 radii
    
    if dset == 2;load('Data/ScratchPeaksXY');   % Scratchbot (ROBIO expts)
        data = ScratchPeaksXY; clear ScratchPeaksXY;% 3 speeds, 3 radii, 8 contacts
    end;
    
    if dset == 3;load('Data/roombaRadiusXY');   % Crunchbot.
        data = roombaRadiusXY; clear roombaRadiusXY;% 4 whiskers, 6 radii, 5 contacts
    end;
    
    % if dset == 4;load('Data/Collated.mat'); % Scratchbot (old expt1, saw/sin, 90-190mm, unpublished)
    % data = expt1.sawradial90deg.xy % DATA NOT READY. Needs re-sorting
    

    smoo = 1; % smooth the data 0|1
    deriv = 1; % take first of second derivatives 0 | 1 | 2
   
    %% PREPROCESS DATA INTO TRAINING + TEST SETS
    
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    % training data, concatenated training data, test data, radius of
    % trial, speed of trial | dataset, label of dataset
    
    
    
    %% Insert loop for reducing training set size
    
  
    
    
    %% %% RUN CLASSIFIERS
    %% Static
    % Specific preprocessing
    
    
    %% Freq
    
    %% Temp
    
    %% Feat
    %  AR_train_all('static')
    
    %% MAP
    % Specific preprocessing
    % Concatenate files
    
    % TRAINING
    
    %   TO DO:  AR_train('SNB'); New universal code?
    %   Old code
    d = linspace(-0.1,0.1,501); d = d(:); % set params for histogram, then rotate
    ns = 10;                              % smooth over 10 samples
    clear logl;
    for c = 1:length(data_train_c);
        logl{c} = AR_train_SNB(data_train_c{c},d,ns);
    end
    
    % TESTING
    
    % Old code
    clear logp;
    for c = 1:length(data_train_c);
        for c1 = 1:length(data_train_c);
            logp{c}(c1) = AR_test_SNB(data_test{c},d,logl{c1});
        end
        [ig,class(c)] = max(squeeze(logp{c}),[],2);
    end

    
    %% GP
    
end % for all datasets