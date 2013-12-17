% This code runs the old classifiers on the 3 different robot datasets,
% under all the relevant conditions for publication in Evans et al 2013
% 'State of the art in whisker-based object localisation'.
% Puts all of the outputs in an array for easy copying into the paper.
% AR_figs.m generates the publcation quality figures. See README.txt for
% more details.

clear all; clc;

% load the four data sets, one at a time.
for dset = 1:3;
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
    % data = expt1.sawradial90deg.xy % DATA NOT READY. Needs re-compiling
    
    % run old scripts to prepare training/test data, editted with counter measures (leave-one-out)
    [data_train,data_test,indRadius,indVelocity] = AR_preprocessing(data,dset)
    
    % Static
    % Specific preprocessing
    
    
    % Freq
    
    % Temp
    
    % Feat
    AR_train_all('static')
    
    % MAP
    % Specific preprocessing
    % Concatenate files
    
    % GP
    
end % for all datasets