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
        curname = dirnames{ndt};
        cd(curname)
        
        clear params model fwhms
        modelfile = dir('*Monotonic-Log-*-Log-*-gFit.mat');
        modelname = modelfile.name;
        newname = [modelname(1:25),'2Lin' ,modelname(26:85), '2Lin', modelname(86:end)]
        load(modelname)
        params.matFileName{1} = [modelfile.folder, '/', newname];
        params.matFileName{2} = newname;
        Log2Lin;
        save(newname, 'params', 'model', 'fwhms');
%         if exist('xvalRefit', 'dir') == 7
%             cd('xvalRefit')
%             clear params model fwhms
%             modelfile = dir('*Monotonic-Log-*-Log-*-fFit.mat');
%             modelname = modelfile.name;
%             newname = [modelname(1:26),'2Lin' ,modelname(27:87), '2Lin', modelname(88:end)]
%             load(modelname)
%             params.matFileName{1} = [modelfile.folder, '/', newname];
%             params.matFileName{2} = newname;
%             Log2Lin;
%             save(newname, 'params', 'model', 'fwhms');
%             cd ../..
%         else
            cd ..
%         end
    end
end