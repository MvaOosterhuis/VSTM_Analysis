%From these ROIs, get data for various models
mapNames=["VMIPS0", "VMIPS1", "VMIPS2", "VMIPS3", "VMIPS4", "VMIPS5", "VMT1", "VMT2", "VMF1", "VMF2", "VMF3"];
% modelProg=["1DGaussian-Lin-", "1DGaussian-Log2Lin-", "Monotonic-Lin", "Monotonic-Log-"]; % "-Lin-",
% modelProg=[ "Monotonic-Lin", "Monotonic-Log-"]; %
% modelProgName = ["LinMonotonic", "LogMonotonic"];%
% modelProgName = ["Log2Lin","LinMonotonic", "LogMonotonic"];%
% whichModel='-gFit-gFit';
% DTs=[13,14,15,16,17];
DTs = [5]; % For Surya Reverse
% DTnames=["Progressive", "Orientation", "Difficulty", "Odd", "Even", "RandomOrder"];
DTnames=["Reverse"];
subs = ["DataS01","DataS13","DataS14","DataS19","DataS16","DataS22"];
ppdirnames = {'/mnt/data/Try_full_run/Working_Memory/S01/FRESH_WorkingMemory_S01_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S13/FRESH_WorkingMemory_S13_All/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S14/FRESH_WorkingMemory_S14/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S19/FRESH_WorkingMemory_S19/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S16/Session1/mrVistaSession',...
    '/mnt/data/Try_full_run/Working_Memory/S22/Session1/mrVistaSession'};


for whichSub = [2]%:length(subs)
    
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
    for whichDT=1%:length(DTs)
        if DTs(whichDT)>0
            VOLUME{1}=viewSet(VOLUME{1}, 'curdt', DTs(whichDT));
            for addModels = 1
                % Find the DTs directory where the model files live
                if addModels == 1
                    folder=[pwd, filesep, 'Gray', filesep, 'Reverse'];
                    modelFile=dir([char(strcat(folder,filesep, '*Log2Lin*-gFit.mat'))]);
                else
                end
                
                
                VOLUME{1}=rmSelect(VOLUME{1}, 1,[folder, filesep, modelFile.name]);
                
                VOLUME{1} = rmLoadDefault(VOLUME{1});
                VOLUME{1}=refreshScreen(VOLUME{1},0);
                
                
                for whichMap=1:length(mapNames)
                    if addModels == 1
                        dataName=char(strcat(subname, '.', mapNames(whichMap), '.Left.Reverse.Log2Lin'));
                    else
                    end
                    
                    
                    %Left hemisphere first entries
                    if okROIs(whichMap,1)
                        dataTmp = RoiDistanceRatio(VOLUME{1}, 2+(whichMap-1)*6, 3+(whichMap-1)*6, 1+(whichMap-1)*6, [], length(VOLUME{1}.ROIs));
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        eval([dataName, '=dataTmp;'])
                    end
                    
                    if addModels == 1
                        dataName=char(strcat(subname, '.', mapNames(whichMap), '.Right.Reverse.Log2Lin'));
                    else
                    end
                    %Right hemisphere
                    if okROIs(whichMap,2)
                        dataTmp = RoiDistanceRatio(VOLUME{1}, 5+(whichMap-1)*6, 6+(whichMap-1)*6, 4+(whichMap-1)*6, [], length(VOLUME{1}.ROIs));
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        eval([dataName, '=dataTmp;'])
                    end
                    
                end
            end
        end
    end
end