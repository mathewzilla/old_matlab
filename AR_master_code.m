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
    if dset == 1;
        fprintf('Loading XY table data\n')
        load('AR_Data/data_XY','data');   % XY table
    end;                                        % 4 sets, 26 speeds, 101 radii
    
    if dset == 2;
        fprintf('Loading Scratchbot data\n')
        load('AR_Data/ScratchPeaksXY');   % Scratchbot (ROBIO expts)
        data = ScratchPeaksXY; clear ScratchPeaksXY;% 3 speeds, 3 radii, 8 contacts
    end;
    
    if dset == 3;
        fprintf('Loading roomba data\n')
        load('AR_Data/roombaRadiusXY');   % Crunchbot.
        data = roombaRadiusXY; clear roombaRadiusXY;% 4 whiskers, 6 radii, 5 contacts
    end;
    
    
    %% PREPROCESS DATA INTO TRAINING + TEST SETS
    smoo = 0; % smooth the data 0|1
    deriv = 0; % take first of second derivatives 0 | 1 | 2 . ATM AR_train_feat.m doesn't work with derivs
    
    fprintf('Preprocessing for template and feature based classifiers... \n')
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    % training data, concatenated training data, test data, radius of
    % trial, speed of trial | dataset, label of dataset
    
    %% %% RUN CLASSIFIERS
    
    %% Temp
    % TRAIN AND TEST CLASSIFIER
    fprintf('Running template based classifier ')
    method = 1;
    clear class;
    [class,t_train{1},t_test{1}] = AR_run_temp(data_train,data_test,indRadius,indVelocity,method);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat{1}] = AR_stats(indRadius,indVelocity,class);
    
    t_train{1}
    t_test{1}
    stat{1}
    
    %% Freq
    fprintf('Running frequency template based classifier ')
    method = 5; % Frequency template
    clear class;
    [class,t_train{2},t_test{2}] = AR_run_temp(data_train,data_test,indRadius,indVelocity,method);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat{2}] = AR_stats(indRadius,indVelocity,class);
    
    t_train{2}
    t_test{2}
    stat{2}
    
    %% Feat
    % TRAIN AND TEST CLASSIFIER
    fprintf('Running feature based classifier ')
    clear class;
    [class,t_train{3},t_test{3}] = AR_run_feat(data_train,data_test,indRadius,indVelocity);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat{3}] = AR_stats(indRadius,indVelocity,class);
    
    t_train{3}
    t_test{3}
    stat{3}
    
    %% MAP
    smoo = 1; deriv = 1;
    fprintf('Preprocessing for NB classifier... \n')
    clear class;
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    
    % TRAIN AND TEST CLASSIFIER
    fprintf('Running Naive Bayes classifier ')
    [class,t_train{4},t_test{4}] = AR_run_SNB(data_train_c,data_test);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat{4}] = AR_stats(indRadius,indVelocity,class);
    
    t_train{4}
    t_test{4}
    stat{4}
    
    %% Static
    fprintf('Running static beam equation based classifier ')
    clear class;
    [class,t_train{5},t_test{5}] = AR_run_static(data_train,data_test,indRadius,indVelocity);
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat{5}] = AR_stats(indRadius,indVelocity,class);
    
    t_train{5}
    t_test{5}
    stat{5}
    
    %% GP
    
    
end % for all datasets

%% TO DO Insert loop for reducing training set size
%% Add 4th dataset
% if dset == 4;load('Data/Collated.mat'); % Scratchbot (old expt1, saw/sin, 90-190mm, unpublished)
% data = expt1.sawradial90deg.xy % DATA NOT READY. Needs re-sorting
