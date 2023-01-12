dirnames = {'ConstantDifficultyAverages (L1)',...
    'ConstantOrientationAverages (L1)',...
    'EvenScansAverages (L1)',...
    'OddScansAverages (L1)',...
    'ProgressiveAverages (L1)',...
    'RandomOrderAverages (L1)'};

ppdirnames = {'/mnt/data/Try_full_run/Working_Memory/S01/FRESH_WorkingMemory_S01_All',...
    '/mnt/data/Try_full_run/Working_Memory/S13/FRESH_WorkingMemory_S13_All',...
    '/mnt/data/Try_full_run/Working_Memory/S14/FRESH_WorkingMemory_S14',...
    '/mnt/data/Try_full_run/Working_Memory/S16/Session1',...
    '/mnt/data/Try_full_run/Working_Memory/S19/FRESH_WorkingMemory_S19'};

for np = 1:length(ppdirnames)
    
    grayname = [ppdirnames{np}, '/mrVistaSession/Gray'];
    cd(grayname)
    
    for ndt = 1:length(dirnames)
        clear params model
        curname = dirnames{ndt};
        cd(curname)
        modelfile = dir('*1DGaussian-FullBlanks*-gFit-gFit.mat');
        modelname = modelfile.name;
        newname = [modelname(1:23),'Lin-' ,modelname(24:82), 'Lin-', modelname(83:end)]
        load(modelname);
        params.matFileName{1} = [modelfile.folder, '/', newname];
        params.matFileName{2} = newname;
        
        save(newname, 'params', 'model');
        cd ..
    end
end