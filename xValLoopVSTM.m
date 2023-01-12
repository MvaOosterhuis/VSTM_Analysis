%%LoopForXVal
close all
clear all

addpath(genpath('/mnt/data/matlabNew/bensfunctions'))
addpath(genpath('/mnt/data/matlabNew/'))
addpath(genpath('/mnt/data/vistasoftOld'))

DTs=[13,14,15,16,17]; %S13, S01, S14, S16, s19

% Define all the models that we can have
models = {
'-Lin-',...
'-Log-'};


ppdirnames = {'/mnt/data/Try_full_run/Working_Memory/S01/FRESH_WorkingMemory_S01_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S13/FRESH_WorkingMemory_S13_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S14/FRESH_WorkingMemory_S14/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S16/Session1/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S19/FRESH_WorkingMemory_S19/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S22/Session1/mrVistaSession'};

paths = ppdirnames;
modelTypes = models;
combinedDT = DTs;

for np = 6%1:length(paths)
    whichSubs = np;    
    % Inside CrossValidate.. combinedDT(4:5) is used for xVal purposes,
    % make sure the order is correct.
    CrossValidateCandidateVSTMModels(paths, whichSubs, combinedDT, modelTypes)
end