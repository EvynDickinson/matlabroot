
% Temp-rate distance tuning curve for a single genotype or across genotypes
% Show how a single genotypes compares in food attraction:temperature for
% different rates of temperature change


%% Select data groups to compare

clear; close all; clc

% Get names of processed data:
baseFolder = getCloudPath;  
structFolder = [baseFolder 'Data structures\'];
list_dirs = dir(structFolder);
list_dirs = {list_dirs(:).name};
list_dirs(1:2) = [];


  

expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
expName = expName(1:end-1);clear expNames




 





%%
% Cross genotype comparisons: 


