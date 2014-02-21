% Code for manipulating AR_Results.mat into the correct format for publication
% Rearranges key results into tables, then prints these tables to file
% To be called from AR_2013 master directory
clear all
load Results/AR_Results.mat

% Option for generating tables, histograms or plots
tables = 0;
hists = 1;

if tables;
    %% DATASET 1 Rearrange for plotting tables
    for i = 1:5;
        results.table1(i,1:4) = results.stat{1,i};
        results.table1(i,5) = results.t_train{1,i};
        results.table1(i,6) = results.t_test{1,i};
    end
    
    % Print results table to .tex file. Code inspired by N.Lawrence
    classy = {'Temp', 'Freq', 'Feat', 'MAP', 'Static'};
    [rows cols] = size(results.table1);
    numSigFigs = 3;
    fid = fopen('Results/results_table1.tex', 'w');
    for i = 1:rows
        fprintf(fid,[classy{i} ' & ']);
        for j = 1:cols
            fprintf(fid, ['$' num2str(results.table1(i, j), numSigFigs) '$']);
            if j < cols
                fprintf(fid, ' & ');
            end
        end
        if i <= rows
            fprintf(fid, '\\\\\n');
        end
    end
    fclose(fid);
    
    %% DATASET 2 Rearrange for plotting tables
    for i = 1:5;
        results.table2(i,1:4) = results.stat{2,i};
        results.table2(i,5) = results.t_train{2,i};
        results.table2(i,6) = results.t_test{2,i};
    end
    
    % Print results table to .tex file.
    
    fid = fopen('Results/results_table2.tex', 'w');
    for i = 1:rows
        fprintf(fid,[classy{i} ' & ']);
        for j = 1:cols
            fprintf(fid, ['$' num2str(results.table2(i, j), numSigFigs) '$']);
            if j < cols
                fprintf(fid, ' & ');
            end
        end
        if i <= rows
            fprintf(fid, '\\\\\n');
        end
    end
    fclose(fid);
    
    %% DATASET 3 Rearrange for plotting tables
    for i = 1:5;
        results.table3(i,1:4) = results.stat{3,i};
        results.table3(i,5) = results.t_train{3,i};
        results.table3(i,6) = results.t_test{3,i};
    end
    
    % Print results table to .tex file.
    [rows cols] = size(results.table3);
    numSigFigs = 3;
    fid = fopen('Results/results_table3.tex', 'w');
    for i = 1:rows
        fprintf(fid,[classy{i} ' & ']);
        for j = 1:cols
            fprintf(fid, ['$' num2str(results.table3(i, j), numSigFigs) '$']);
            if j < cols
                fprintf(fid, ' & ');
            end
        end
        if i <= rows
            fprintf(fid, '\\\\\n');
        end
    end
    fclose(fid);
end


%% Plotting histograms of errors
if hists;
    k = 0;
    for i = 1:3;
        for j = 1:5;
            k = k + 1;
            figure(1);
            subplot(3,5,k);
            m1 = results.errorR{i,j};
            hist(m1)
            axis square
            
            figure(2);
            subplot(3,5,k);
            m1 = results.errorR{i,j};
            hist(m1)
            axis square
            
        end
    end
end