%From these ROIs, get data for various models
mapNames=["VMIPS0", "VMIPS1", "VMIPS2", "VMIPS3", "VMIPS4", "VMIPS5", "VMT1", "VMT2", "VMF1", "VMF2", "VMF3"];
% modelProg=["1DGaussian-Log2Lin-"]; % "-Lin-",
modelProg=["1DGaussian-Lin-", "1DGaussian-Log2Lin-", "Monotonic-Lin", "Monotonic-Log-"]; % "-Lin-",
% modelProg=[ "Monotonic-Lin", "Monotonic-Log-"]; %
% modelProgName = ["LinMonotonic", "LogMonotonic"];%
% modelProgName = ["Log2Lin"];%
modelProgName = ["Lin","Log2Lin","LinMonotonic", "LogMonotonic"];%
% whichModel='-gFit-gFit';
DTs=[13,14,15,16,17,18];
DTnames=["Progressive", "Orientation", "Difficulty", "Odd", "Even", "RandomOrder"];
subs = ["DataS01","DataS13","DataS14","DataS19","DataS16","DataS22"];
ppdirnames = {'/mnt/data/Try_full_run/Working_Memory/S01/FRESH_WorkingMemory_S01_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S13/FRESH_WorkingMemory_S13_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S14/FRESH_WorkingMemory_S14/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S19/FRESH_WorkingMemory_S19/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S16/Session1/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S22/Session1/mrVistaSession'};


for whichSub = [1:6]%:length(subs)
    
    if whichSub == 5
        DTs = DTs(1:5);
        DTnames = DTnames(1:5);
    end
    
    if whichSub == 4
        modelProg=["1DGaussian-Lin-", "1DGaussian-Log2Lin-"]; % "-Lin-",
        modelProgName = ["Lin","Log2Lin"];%
    else
        modelProg=["1DGaussian-Lin-", "1DGaussian-Log2Lin-", "Monotonic-Lin", "Monotonic-Log-"]; % "-Lin-",
        modelProgName = ["Lin","Log2Lin","LinMonotonic", "LogMonotonic"];%
    end
    
    subdir = ppdirnames{whichSub};
    subname = subs(whichSub);
    cd(subdir)
    
    close all
    VOLUME{1} = mrVista('3')
    VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
    okROIs=zeros(length(mapNames), 2);
    for whichMap=1:length(mapNames)
        [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'Map')), [],[],[],1);
        if ok
            [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'Low')), [],[],[],1);
            [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'High')), [],[],[],1);
            okROIs(whichMap,1)=ok;
        else
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 1]);
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 2]);
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 3]);
        end
        
        [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'Map')), [],[],[],1);
        if ok
            [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'Low')), [],[],[],1);
            [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'High')), [],[],[],1);
            okROIs(whichMap,2)=ok;
        else
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 4]);
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 5]);
            VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 6]);
        end
    end
    VOLUME{1} = loadROI(VOLUME{1}, 'gray-Layer1', [],[],[],1)
    VOLUME{1} = refreshScreen(VOLUME{1}, 0);
    for whichDT=1:length(DTs)
        if DTs(whichDT)>0
            VOLUME{1}=viewSet(VOLUME{1}, 'curdt', DTs(whichDT));
            for whichProg=1:length(modelProg)
                % Find the DTs directory where the model files live
                folder=[pwd filesep 'Gray' filesep dataTYPES(DTs(whichDT)).name];
                
                modelFile=dir([char(strcat(folder,filesep, '*', modelProg(whichProg), '*-gFit.mat'))]);
                VOLUME{1}=rmSelect(VOLUME{1}, 1,[folder, filesep, modelFile.name]);
                
                VOLUME{1} = rmLoadDefault(VOLUME{1});
                VOLUME{1}=refreshScreen(VOLUME{1},0);
                
                if (whichDT==4 || whichDT==5) && whichProg < 3
                    xvalFolder=[pwd filesep 'Gray' filesep dataTYPES(DTs(9-whichDT)).name filesep '/xvalRefit'];
                    xvalModelFile=dir([char(strcat(xvalFolder,filesep, '*', modelProg(whichProg), '*-fFit.mat'))]);
                    xvalModel=load([xvalFolder, filesep, xvalModelFile.name], 'model');
                    xvalVes=rmGet(xvalModel.model{1}, 've');
                elseif (whichDT==2 || whichDT==3) && whichProg < 3
                    xvalFolder=[pwd filesep 'Gray' filesep dataTYPES(DTs(5-whichDT)).name filesep '/xvalRefit'];
                    xvalModelFile=dir([char(strcat(xvalFolder,filesep, '*', modelProg(whichProg), '*-fFit.mat'))]);
                    xvalModel=load([xvalFolder, filesep, xvalModelFile.name], 'model');
                    xvalVes=rmGet(xvalModel.model{1}, 've');
                end
                
                for whichMap=1:length(mapNames)
                    dataName=char(strcat(subname, '.', mapNames(whichMap), '.Left.', modelProgName(whichProg), '.', DTnames(whichDT)));
                    
                    %Left hemisphere first entries
                    if whichMap == 4
                        breakcheck = 1;
                    end
                    if okROIs(whichMap,1)
                        dataTmp = RoiDistanceRatio(VOLUME{1}, 2+(whichMap-1)*6, 3+(whichMap-1)*6, 1+(whichMap-1)*6, [], length(VOLUME{1}.ROIs));
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        
                        if (whichDT==2 || whichDT==3 || whichDT==4 || whichDT==5) && whichProg < 3
                            dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                            dataTmp.folderXval=xvalFolder;
                            dataTmp.modelFileXval=xvalModelFile.name;
                        end
                        
                        eval([dataName, '=dataTmp;'])
                    end
                    
                    dataName=char(strcat(subname, '.', mapNames(whichMap), '.Right.', modelProgName(whichProg), '.', DTnames(whichDT)));
                    %Right hemisphere
                    if okROIs(whichMap,2)
                        dataTmp = RoiDistanceRatio(VOLUME{1}, 5+(whichMap-1)*6, 6+(whichMap-1)*6, 4+(whichMap-1)*6, [], length(VOLUME{1}.ROIs));
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        if (whichDT==2 || whichDT==3 || whichDT==4 || whichDT==5) && whichProg < 3
                            dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                            dataTmp.folderXval=xvalFolder;
                            dataTmp.modelFileXval=xvalModelFile.name;
                        end
                        eval([dataName, '=dataTmp;'])
                    end
                    
                end
            end
        end
    end
end