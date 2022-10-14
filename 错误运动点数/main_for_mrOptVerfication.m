%% 验证优化方案
clear;
clc;
close all;
%%
global count_sub; % 被试计数
count_sub = 0;
data_foldername = 'data/';
functions_foldername = 'functions/';
addpath(functions_foldername);
subdir = dir(data_foldername);
for idxDir = 1:length(subdir)
    if(isequal(subdir(idxDir).name, '.') || ...
            isequal(subdir(idxDir).name, '..'))
        continue;
    end
    if subdir(idxDir).isdir
        GUIdatapath = fullfile(data_foldername,subdir(idxDir).name);
        dat = dir(GUIdatapath);
        for i = 1:length(dat)
            filenameTemp = strsplit(dat(i).name, '.');
            if isequal(filenameTemp{end}, 'mat')
                dataname = strcat(filenameTemp{1}, '.mat');
                fileLoadPath = fullfile(GUIdatapath, dataname);
                GUIdataTemp = load(fileLoadPath);
                %%
                field = fieldnames(GUIdataTemp);
                GUIdata = GUIdataTemp.(field{1});
                algorithmForWaveOpt(GUIdata,fileLoadPath,filenameTemp{1}); %
            end
        end
    else
        filenameTemp = strsplit(subdir(idxDir).name, '.');
        if isequal(filenameTemp{end}, 'mat')
            fileLoadPath = fullfile(data_foldername,subdir(idxDir).name);
            GUIdataTemp = load(fileLoadPath);
            field = fieldnames(GUIdataTemp);
            GUIdata = GUIdataTemp.(field{1});
            count_sub = count_sub + 1;
            algorithmForWaveOpt(GUIdata,fileLoadPath,filenameTemp{1}); %
            %             try
            %                 algorithmForWaveOpt(GUIdata,fileLoadPath,filenameTemp{1}); %
            %             catch exception
            %                 continue;
            %             end
        end
    end
end