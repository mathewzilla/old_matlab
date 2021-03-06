% This code runs the old classifiers on the 3 different robot datasets,
% under all the relevant conditions for publication in Evans et al 2014
% 'State of the art in whisker-based object localisation'.
% Puts all of the outputs in an array for easy copying into the paper.
% AR_figs.m generates the publcation quality figures. See README.txt for
% more details.

% clear all; clc;

% load the three (not four) data sets, one at a time.
for dset = 1:3;% :3; other numbers don't work yet
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
        % NOTE: no velocity differences, but 4 separate whiskers.
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
    clear class t_train t_test;
    [class,t_train,t_test] = AR_run_temp(data_train,data_test,indRadius,indVelocity,method);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class);
    results.outR{dset,1} = outR;
    results.outV{dset,1} = outV;
    results.errorR{dset,1} = errorR;
    results.errorV{dset,1} = errorV;
    results.t_train{dset,1} = t_train;
    results.t_test{dset,1} = t_test;
    results.stat{dset,1} = stat;
    clear outR outV errorR errorV stat;
    
    %% Freq
    fprintf('Running frequency template based classifier ')
    method = 5; % Frequency template
    clear class t_train t_test;
    [class,t_train,t_test] = AR_run_temp(data_train,data_test,indRadius,indVelocity,method);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class);
    
    results.outR{dset,2} = outR;
    results.outV{dset,2} = outV;
    results.errorR{dset,2} = errorR;
    results.errorV{dset,2} = errorV;
    results.t_train{dset,2} = t_train;
    results.t_test{dset,2} = t_test;
    results.stat{dset,2} = stat;
    clear outR outV errorR errorV stat;
    
    %% Feat
    % TRAIN AND TEST CLASSIFIER
    fprintf('Running feature based classifier ')
    clear class t_train t_test;
    [class,t_train,t_test] = AR_run_feat(data_train,data_test,indRadius,indVelocity,dset);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class);
    
    results.outR{dset,3} = outR;
    results.outV{dset,3} = outV;
    results.errorR{dset,3} = errorR;
    results.errorV{dset,3} = errorV;
    results.t_train{dset,3} = t_train;
    results.t_test{dset,3} = t_test;
    results.stat{dset,3} = stat;
    clear outR outV errorR errorV stat;
    
    %% MAP
    smoo = 1; deriv = 1;
    fprintf('Preprocessing for NB classifier... \n')
    clear class t_train t_test;
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    
    % TRAIN AND TEST CLASSIFIER
    fprintf('Running Naive Bayes classifier ')
    [class,t_train,t_test] = AR_run_SNB(data_train_c,data_test);
    
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    [outR,outV,errorR,errorV,stat] = AR_stats(indRadius,indVelocity,class);
    
    results.outR{dset,4} = outR;
    results.outV{dset,4} = outV;
    results.errorR{dset,4} = errorR;
    results.errorV{dset,4} = errorV;
    results.t_train{dset,4} = t_train;
    results.t_test{dset,4} = t_test;
    results.stat{dset,4} = stat;
    clear outR outV errorR errorV stat;
    
    %% Static
    smoo = 0; deriv = 0;
    fprintf('Preprocessing for static classifier... \n')
    clear class t_train t_test;
    [data_train,data_train_c,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset,smoo,deriv);
    
    fprintf('Running static beam equation based classifier ')
    clear class;
    [class,t_train,t_test] = AR_run_static(data_train,data_test,indRadius,indVelocity,dset);
    % COMPUTE RESULTS
    fprintf('Computing results \n')
    % STAT NOT WORKING FOR STATIC AS CLASSIFIER ISN'T ACCURATE
    % [outR,outV,errorR,errorV,stat{dset,5}] = AR_stats(indRadius,indVelocity,class);
    outR = class(1,:);
    outV = class(2,:);
    errorR = outR - indRadius(1,:);
    errorV = outV - indVelocity(1,:);
    
    stat(1) = mean(abs(errorV)); % TO DO: Multiply by actual speed units
    stat(2) = std(abs(errorV));
    stat(3) = mean(abs(errorR));
    stat(4) = std(abs(errorR));
    
    results.outR{dset,5} = outR;
    results.outV{dset,5} = outV;
    results.errorR{dset,5} = errorR;
    results.errorV{dset,5} = errorV;
    results.t_train{dset,5} = t_train;
    results.t_test{dset,5} = t_test;
    results.stat{dset,5} = stat;
    clear outR outV errorR errorV stat;
    
    %% GP
    
    
end % for all datasets

% savefile = 'Results/AR_Results.mat';
% save(savefile,'results')

%% TO DO Insert loop for reducing training set size
%% Add 4th dataset
% if dset == 4;load('Data/Collated.mat'); % Scratchbot (old expt1, saw/sin, 90-190mm, unpublished)
% data = expt1.sawradial90deg.xy % DATA NOT READY. Needs re-sorting
