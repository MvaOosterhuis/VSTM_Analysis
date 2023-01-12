% Avoid confusing poor mrVista/MATlab

addpath(genpath('/mnt/data/matlabNew/vistasoft_preproc')) % For roiRestrictByLayer

close all; 

VOLUME{1} = mrVista('3');
fprintf('Reducing tSeries for dataTYPE TimingSweeps to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ProgressiveAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE FirstHalf to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ConstantOrientationAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE SecondHalf to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','ConstantDifficultyAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE OddScans to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','OddScansAverages');
collapseOverLayersTseries(VOLUME{1});

fprintf('Reducing tSeries for dataTYPE EvenScans to Layer-1. \n')
VOLUME{1}=viewSet(VOLUME{1},'curdt','EvenScansAverages');
collapseOverLayersTseries(VOLUME{1});


fprintf('Making and Saving the Gray ROI. \n.')
%Create an ROI of the layer 1 voxels and save to local file
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
% VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI);
VOLUME{1} = roiRestrictByLayer(VOLUME{1},VOLUME{1}.selectedROI,1);
VOLUME{1} = refreshScreen(VOLUME{1},0);
saveROI(VOLUME{1}, 'selected', 1);