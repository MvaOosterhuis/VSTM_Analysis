dirnames = {'ConstantDifficultyAverages (L1)',...
    'ConstantOrientationAverages (L1)',...
    'EvenScansAverages (L1)',...
    'OddScansAverages (L1)',...
    'ProgressiveAverages (L1)',...
    'RandomOrderAverages (L1)'};

ppdirnames = {'/mnt/data/Try_full_run/Working_Memory/S01/FRESH_WorkingMemory_S01_All',...
    '/mnt/data/Try_full_run/Working_Memory/S13/FRESH_WorkingMemory_S13_All',...
    '/mnt/data/Try_full_run/Working_Memory/S14/FRESH_WorkingMemory_S14',...
    '/mnt/data/Try_full_run/Working_Memory/S19/FRESH_WorkingMemory_S19',...
    '/mnt/data/Try_full_run/Working_Memory/S16/Session1'};

for np = 1:length(ppdirnames)
    if np == 5
        dirnames = dirnames(1:5);
    end
    
    grayname = [ppdirnames{np}, '/mrVistaSession/Gray'];
    cd(grayname)
    
    for ndt = 1:length(dirnames)
        clear params model
        curname = dirnames{ndt};
        cd(curname)
        modelfile = dir('*Monotonic-FullBlanks*-gFit.mat');
        modelname = modelfile.name;
        newname = [modelname(1:22),'Lin-' ,modelname(23:80), 'Lin-', modelname(81:end)]
        load(modelname);
        params.matFileName{1} = [modelfile.folder, '/', newname];
        params.matFileName{2} = newname;
        
        save(newname, 'params', 'model');
        cd ..
    end
end