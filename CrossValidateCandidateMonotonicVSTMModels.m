function CrossValidateCandidateMonotonicVSTMModels(paths, whichSubs, combinedDT, modelTypes)%, ParamsPathMix, ParamsPathSessions)
%Cross validates a set of timing selectivity pRF models in a specified
%folder

%BMH, 01/2020

for thisSub=whichSubs
    cd(paths{thisSub})
    VOLUME{1}=mrVista('3');
    allXvalDTs=combinedDT(4:5);
    load('mrSESSION.mat', 'dataTYPES')
    
    % First make sure the correct dirs are selected, and loop over ALL of
    % them
    for nm = 1:length(modelTypes)
        modelType = modelTypes{nm};
        
        % Then start even/odd stuff
        for n=1:length(allXvalDTs)
            files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/*', modelType ,'*-gFit.mat']);
            thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/'];
            otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/'];
            if ~exist([otherPath,'/xval'], 'dir')
                eval(['!mkdir ',  '"',otherPath, 'xval"']);
                eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
            else
            end
            
            for whichFile=1:length(files)
                eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
            end
        end
        
        % Loop for xVal
        for ldt = 1:length(allXvalDTs)
            dtname{ldt} = dataTYPES(allXvalDTs(ldt)).name;
        end
        for whichDT=1:length(allXvalDTs)
            origName = [pwd '/Gray/' dtname{whichDT}];
            origmodelFile = dir([origName '/*', modelType, '*-gFit.mat']);
            folderName=[pwd '/Gray/' dtname{whichDT}, '/xval'];
            refitName=[pwd '/Gray/' dtname{whichDT}, '/xvalRefit'];
            modelFiles=dir([folderName '/*', modelType, '*-gFit.mat']);
            for whichModel=1:length(modelFiles)
                % Actual xVal call
                try
                    load([origName, '/', origmodelFile.name])
                    origmodel = model;
                    load([folderName, '/', modelFiles.name])
                    xvalmodel = model;
                    
                    % Define indices
                    xvals = xvalmodel{1}.x0;
                    origvals = origmodel{1}.x0;
                    newvals = xvals.*origvals; %(anything*0 is 0, both 1 == 1, both 2 == 4, any difference == 2)
                    inds = newvals == 1 | newvals == 4;
                    
                    %Use indices in xval model to remake it
                    xvalmodel{1}.rawrss(~inds) = 0;
                    xvalmodel{1}.rss(~inds) = 0;
                    xvalmodel{1}.x0(~inds) = 0;
                    
                    model = xvalmodel;
                    refitfilename = [modelFiles.name(1:end-4), '-fFit.mat'];
                    save([refitName, '/', refitfilename],'model', 'params')
                catch
                    fprintf([
                        '#######################################################################\n',...
                        '#######################################################################\n',...
                        'Error during PostSearch in model: %s.\n',...
                        '#######################################################################\n',...
                        '#######################################################################\n']);
                    continue;
                end
            end
        end
    end
    close(1)
    
    % Also cross validate Orientation with Difficulty, where we expect the
    % two to be relatively uncorrelated if they indeed lead to different
    % response. A good correlation means both good model and good overlap
    % between the task conditions
    VOLUME{1}=mrVista('3');
    allXvalDTs=combinedDT(2:3);
    load('mrSESSION.mat', 'dataTYPES')
    
    % First make sure the correct dirs are selected, and loop over ALL of
    % them
    for nm = 1:length(modelTypes)
        modelType = modelTypes{nm};
        
        % Then start even/odd stuff
        for n=1:length(allXvalDTs)
            files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/*', modelType ,'*-gFit.mat']);
            thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/'];
            otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/'];
            if ~exist([otherPath,'/xval'], 'dir')
                eval(['!mkdir ',  '"',otherPath, 'xval"']);
                eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
            else
            end
            
            for whichFile=1:length(files)
                eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
            end
        end
        
        % Loop for xVal
        for ldt = 1:length(allXvalDTs)
            dtname{ldt} = dataTYPES(allXvalDTs(ldt)).name;
        end
        for whichDT=1:length(allXvalDTs)
            origName = [pwd '/Gray/' dtname{whichDT}];
            origmodelFile = dir([origName '/*', modelType, '*-gFit.mat']);
            folderName=[pwd '/Gray/' dtname{whichDT}, '/xval'];
            refitName=[pwd '/Gray/' dtname{whichDT}, '/xvalRefit'];
            modelFiles=dir([folderName '/*', modelType, '*-gFit.mat']);
            for whichModel=1:length(modelFiles)
                % Actual xVal call
                try
                    load([origName, '/', origmodelFile.name])
                    origmodel = model;
                    load([folderName, '/', modelFiles.name])
                    xvalmodel = model;
                    
                    % Define indices
                    xvals = xvalmodel{1}.x0;
                    origvals = origmodel{1}.x0;
                    newvals = xvals.*origvals; %(anything*0 is 0, both 1 == 1, both 2 == 4, any difference == 2)
                    inds = newvals == 1 | newvals == 4;
                    
                    %Use indices in xval model to remake it
                    xvalmodel{1}.rawrss(~inds) = 0;
                    xvalmodel{1}.rss(~inds) = 0;
                    xvalmodel{1}.x0(~inds) = 0;
                    
                    model = xvalmodel;
                    refitfilename = [modelFiles.name(1:end-4), '-fFit.mat'];
                    save([refitName, '/', refitfilename],'model', 'params')
                catch
                    fprintf([
                        '#######################################################################\n',...
                        '#######################################################################\n',...
                        'Error during PostSearch in model: %s.\n',...
                        '#######################################################################\n',...
                        '#######################################################################\n']);
                    continue;
                end
            end
        end
    end
    close(1)
end
