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
        fprintf('Loading XY table data\n')
    end;                                        % 4 sets, 26 speeds, 101 radii
    
    if dset == 2;load('Data/ScratchPeaksXY');   % Scratchbot (ROBIO expts)
        data = ScratchPeaksXY; clear ScratchPeaksXY;% 3 speeds, 3 radii, 8 contacts
        fprintf('Loading Scratchbot data\n')
    end;
    
    if dset == 3;load('Data/roombaRadiusXY');   % Crunchbot.
        data = roombaRadiusXY; clear roombaRadiusXY;% 4 whiskers, 6 radii, 5 contacts
        fprintf('Loading roomba data\n')
    end;

    smoo = 1; % smooth the data 0|1
    deriv = 1; % take first of second derivatives 0 | 1 | 2
   
    %% PREPROCESS DATA INTO TRAINING + TEST SETS
    
     fprintf('Preprocessing... \n')
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    % training data, concatenated training data, test data, radius of
    % trial, speed of trial | dataset, label of dataset
   
    %% %% RUN CLASSIFIERS
    
    %% Static
    
    %% Freq
    
    %% Temp
    
    %% Feat
    
    %% MAP
    % TRAIN AND TEST CLASSIFIER
     fprintf('Running Naive Bayes classifier ')
    [class,t_train,t_test] = AR_run_SNB(data_train_c,data_test);

    % COMPUTE RESULTS
     fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class);

    t_train
    t_test
    stat
    %% GP
    
end % for all datasets

    %% TO DO Insert loop for reducing training set size
    %% Add 4th dataset    
    % if dset == 4;load('Data/Collated.mat'); % Scratchbot (old expt1, saw/sin, 90-190mm, unpublished)
    % data = expt1.sawradial90deg.xy % DATA NOT READY. Needs re-sorting
    