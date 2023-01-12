% DTs=[13,14,15,16,17]; %S13, S01, S16
% % DTs=[15,16,17,18,19]; % S14, S19
%S13 and S01, DTs: 13 = progressive; 14 cOrientation; 15 cDifficulty; 16&17 Odd and Even; 
% 18 = RandomOrder
DTs = [23];

% Make sure version of vistasoft is correct
oldversion= '/mnt/data/matlabNew/vistasoft_preproc';
newversion= '/mnt/data/vistasoftOld';
rmpath(genpath(oldversion))
addpath(genpath(newversion))

% Set the right directory
ppdir = pwd;

% Load correct parameters for modeling
load('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat')
% Progressive Log
logVSTMProgressive.paramsFile= [ppdir, '/Stimuli/params_log_WM_15_progressive.mat'];
logVSTMProgressive.imFile= [ppdir, '/Stimuli/None'];
logVSTMProgressive.jitterFile= [ppdir, '/Stimuli/None'];
% Progressive Linear
linVSTMProgressive.paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_progressive.mat'];
linVSTMProgressive.imFile= [ppdir, '/Stimuli/None'];
linVSTMProgressive.jitterFile= [ppdir, '/Stimuli/None'];

% Select relevant ROI
roi{1} = 'gray-Layer1';

% Open session
VOLUME{1} = mrVista('3');

% Run Progressive Log
setAllRetParams(logVSTMProgressive, DTs);
rmRunVSTMScriptLog(VOLUME{1},DTs, roi{1},1)
close all

% Open Session
VOLUME{1} = mrVista('3');

% Run Progressive Linear
setAllRetParams(linVSTMProgressive, DTs);
rmRunVSTMScriptLinear(VOLUME{1},DTs, roi{1},1)
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RandomOrder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change DTs
% DTs = 18; % s13, s01, s16, s14
DTs = 19; % S22

% Reload correct parameters
load('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat')
ppdir = pwd;

% Init
for n = 1:6
    linVSTMRandom(n).paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_', num2str(n) ,'random.mat'];
    linVSTMRandom(n).imFile= [ppdir, '/Stimuli/None'];
    linVSTMRandom(n).jitterFile= [ppdir, '/Stimuli/None'];
    logVSTMRandom(n).paramsFile= [ppdir, '/Stimuli/params_log_WM_15_', num2str(n) ,'random.mat'];
    logVSTMRandom(n).imFile= [ppdir, '/Stimuli/None'];
    logVSTMRandom(n).jitterFile= [ppdir, '/Stimuli/None'];
end

% Random Linear
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(linVSTMRandom, DTs);
% rmRunVSTMScriptLinear(VOLUME{1}, DTs, roi{1}, wSearch) % wSearch
rmRunVSTMScriptLinear(VOLUME{1}, DTs, 'gray-Layer1', wSearch) % wSearch
close all

% Random Log
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(logVSTMRandom, DTs);
rmRunVSTMScriptLog(VOLUME{1}, DTs, 'gray-Layer1', wSearch)
close all
%%
% Odd
linVSTMOdd = linVSTMRandom([1,3,5]);
DTs = [21];
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(linVSTMOdd, DTs);
rmRunVSTMScriptLinear(VOLUME{1}, DTs, 'gray-Layer1', wSearch) % wSearch
close all

logVSTMOdd = logVSTMRandom([1,3,5]);
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(logVSTMOdd, DTs);
rmRunVSTMScriptLog(VOLUME{1}, DTs, 'gray-Layer1', wSearch)
close all

% Even
linVSTMEven = linVSTMRandom([2,4,6]);
DTs = [22];
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(linVSTMEven, DTs);
rmRunVSTMScriptLinear(VOLUME{1}, DTs, 'gray-Layer1', wSearch) % wSearch
close all

logVSTMEven = logVSTMRandom([2,4,6]);
VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf.
setAllRetParams(logVSTMEven, DTs);
rmRunVSTMScriptLog(VOLUME{1}, DTs, 'gray-Layer1', wSearch)
close all



%% For Combined Session
% Combined Log (1:6 random, 7 progressive)
DTs = 4;

load('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat')
ppdir = pwd;
logVSTMBoth = logVSTMRandom;
logVSTMBoth(7) = logVSTMProgressive;
for n = 1:6
    logVSTMBoth(n).paramsFile= [ppdir, '/Stimuli/params_log_WM_15_', num2str(n) ,'random.mat'];
    logVSTMBoth(n).imFile= [ppdir, '/Stimuli/None'];
    logVSTMBoth(n).jitterFile= [ppdir, '/Stimuli/None'];
end
logVSTMBoth(7).paramsFile= [ppdir, '/Stimuli/params_log_WM_15_progressive.mat'];
logVSTMBoth(7).imFile= [ppdir, '/Stimuli/None'];
logVSTMBoth(7).jitterFile= [ppdir, '/Stimuli/None'];

VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf. Hrf still breaks though
setAllRetParams(logVSTMBoth, DTs);
rmRunVSTMScriptLog(VOLUME{1}, DTs, 'gray-Layer1', wSearch)

%% Random Odd/Even

% randomOdd = [1,3,5,6,4,2,3,2,1]; %S1
% randomEven = [2,4,6,5,3,1,4,6,5]; %S1


%%
% Lin
% Combined Log (1:6 random, 7 progressive)
DTs = 4;

load('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat')
ppdir = pwd;
linVSTMBoth = linVSTMRandom;
linVSTMBoth(7) = linVSTMProgressive;
for n = 1:6
    linVSTMBoth(n).paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_', num2str(n) ,'random.mat'];
    linVSTMBoth(n).imFile= [ppdir, '/Stimuli/None'];
    linVSTMBoth(n).jitterFile= [ppdir, '/Stimuli/None'];
end
linVSTMBoth(7).paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_progressive.mat'];
linVSTMBoth(7).imFile= [ppdir, '/Stimuli/None'];
linVSTMBoth(7).jitterFile= [ppdir, '/Stimuli/None'];

VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf. Hrf still breaks though
setAllRetParams(linVSTMBoth, DTs);
rmRunVSTMScriptLinear(VOLUME{1}, DTs, 'gray-Layer1', wSearch)

%% For the AllScans condition
% S13: 1:16 is progressive, 17:34 is random
DTs = 20;

load('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat')
ppdir = pwd;

% Init
linVSTMBoth(1:16) = linVSTMProgressive;
linVSTMBoth(17:34) = linVSTMRandom(1);
logVSTMBoth(1:16) = logVSTMProgressive;
logVSTMBoth(17:34) = logVSTMRandom(1);

randomordervec = [6,1,3,2,4,5,4,5,2,3,6,1,2,4,1,5,6,3]; %S13
% randomordervec = [1,2,3,4,5,6,6,5,4,3,2,1,3,4,2,6,1,5]; %S1

% Specify Progressive
for n = 1:16
    linVSTMBoth(n).paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_progressive.mat'];
    linVSTMBoth(n).imFile= [ppdir, '/Stimuli/None'];
    linVSTMBoth(n).jitterFile= [ppdir, '/Stimuli/None'];
end

% Specify Random
for n = 1:18
    linVSTMBoth(n+16).paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_', num2str(randomordervec(n)) ,'random.mat'];
    linVSTMBoth(n+16).imFile= [ppdir, '/Stimuli/None'];
    linVSTMBoth(n+16).jitterFile= [ppdir, '/Stimuli/None'];
end

% Specify Progressive
for n = 1:16
    logVSTMBoth(n).paramsFile= [ppdir, '/Stimuli/params_log_WM_15_progressive.mat'];
    logVSTMBoth(n).imFile= [ppdir, '/Stimuli/None'];
    logVSTMBoth(n).jitterFile= [ppdir, '/Stimuli/None'];
end

% Specify Random
for n = 1:18
    logVSTMBoth(n+16).paramsFile= [ppdir, '/Stimuli/params_log_WM_15_', num2str(randomordervec(n)) ,'random.mat'];
    logVSTMBoth(n+16).imFile= [ppdir, '/Stimuli/None'];
    logVSTMBoth(n+16).jitterFile= [ppdir, '/Stimuli/None'];
end


VOLUME{1} = mrVista('3');
wSearch = 1; % 1 for single grid fit, 8 for grid, hrf and grid using hrf. Hrf still breaks though

modelTypes = [1:2];
% ncores=length(modelTypes);
% maxCores = 4;
% if ncores>maxCores
%     ncores=maxCores;
% end

for LinLog = modelTypes
    if LinLog == 1 % Linear case
        
        setAllRetParams(linVSTMBoth, DTs);
        rmRunVSTMScriptLinear(VOLUME{1}, DTs, 'gray-Layer1', wSearch)
        
    else % log case
        
        setAllRetParams(logVSTMBoth, DTs);
        rmRunVSTMScriptLog(VOLUME{1}, DTs, 'gray-Layer1', wSearch)
    end
end
% delete(gcp('nocreate'))

%% Make stimuli for random conditions

ppdir = pwd;

linRandom = linVSTMProgressive;
linRandom.nCycles = 3;
linRandom.nUniqueRep = 3;
linRandom.prescanDuration = 9;
linRandom.nDCT = 0.5;
linRandom.nFrames = 324;


for n = 1:6
    
    linRandom.paramsFile= [ppdir, '/Stimuli/params_lin_WM_15_', num2str(n) ,'random.mat'];
    linRandom.imFile= [ppdir, '/Stimuli/None'];
    linRandom.jitterFile= [ppdir, '/Stimuli/None'];
    
    linVSTMRandom(n) = linRandom;
end


% And Log

linRandom = logVSTMProgressive;
linRandom.nCycles = 3;
linRandom.nUniqueRep = 3;
linRandom.prescanDuration = 9;
linRandom.nDCT = 0.5;
linRandom.nFrames = 324;

for n = 1:6
    
    linRandom.paramsFile= [ppdir, '/Stimuli/params_log_WM_15_', num2str(n) ,'random.mat'];
    linRandom.imFile= [ppdir, '/Stimuli/None'];
    linRandom.jitterFile= [ppdir, '/Stimuli/None'];
    
    logVSTMRandom(n) = linRandom;
end

save('/mnt/data/Try_full_run/Working_Memory/VSTMparams.mat', 'linVSTMProgressive', 'linVSTMRandom', 'logVSTMProgressive', 'logVSTMRandom')

%% Before running

% dtnames = {'','','','',''}
% 
% VOLUME{1}=viewSet(VOLUME{1},'curdt','SecondHalf');
% collapseOverLayersTseries(VOLUME{1});

close all; 

VOLUME{1} = mrVista('3');
fprintf('Reducing tSeries for dataTYPE ProgressiveAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ProgressiveAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE ConstantOrientationAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ConstantOrientationAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE ConstantDifficultyAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ConstantDifficultyAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE ConstantDifficultyAverages to Layer-1. \n')
VOLUME{2}=viewSet(VOLUME{2},'curdt','SmallDifficultyAverages');
collapseOverLayersTseries(VOLUME{2});

fprintf('Reducing tSeries for dataTYPE OddScansAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','OddScansAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE EvenScansAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','EvenScansAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE RandomOrderAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','RandomOrderAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE RandomOrderAverages to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ReverseAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Done collapsing over layers. \n')

%% If gray isn't made yet
VOLUME{1} = mrVista('3');
fprintf('Making and Saving the Gray ROI. \n.')
%Create an ROI of the layer 1 voxels and save to local file
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
% VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI);
VOLUME{1} = roiRestrictByLayer(VOLUME{1},VOLUME{1}.selectedROI,1);
VOLUME{1} = refreshScreen(VOLUME{1},0);
saveROI(VOLUME{1}, 'selected', 1);
