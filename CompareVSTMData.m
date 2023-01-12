% CompareVSTMData

%Plots bar chart and pairwise differences between candidate model fits.
clear barData targetDataOdd targetDataEven

bfp = "/mnt/data/matlab/bensfunctions_complete/";
vsp = "/mnt/data/vistasoftOld";
sp = "/mnt/data/";
sp2 = "/mnt/data/Try_full_run/";

addpath(genpath(bfp));
addpath(genpath(vsp));
addpath(sp);
addpath(sp2);

loadallnames;
% mapNames = ["VMIPS0","VMIPS1","VMIPS2","VMIPS3","VMIPS4","VMIPS5","VMT1","VMT2","VMF1","VMF2","VMF3"];
mapNames = ["VMIPS0","VMT1","VMT2","VMIPS1","VMIPS2","VMIPS3","VMIPS4","VMIPS5","VMF3","VMF2","VMF1"];
% subjectOrder = {'DataS01','DataS13','DataS14','DataS16','DataS19'};
subjectOrder = {'DataS01','DataS13','DataS14','DataS16', 'DataS19', 'DataS22'};
hemispheres = {'Left', 'Right'};
modelNamesAll =     {'Lin','Log2Lin', 'LinMonotonic', 'LogMonotonic'};
% modelNamesAll = {'Lin','Log2Lin'};

%use 1-ratio!!
%% Model Comparison WITHIN ROI
thr= 0.1;

%
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix
%Determine where to find data in structure
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1:length(modelNamesAll)
                % for variance explained
                targetDataOdd{whichSub, whichMap, whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Odd.vesXval'));
                targetDataEven{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Even.vesXval'));
            end
        end
    end
end

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1:length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataOdd{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataEven{whichSub, whichMap, whichHemi, modelN});
                    end
                    for oddEven=1:2
                        indices=max(barData{whichSub, whichMap, whichHemi, oddEven},[],2)>=0.2;
                        barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=median(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
                    end
                end
            end
        end
    end
end

%Compute and plot means and confidence intervals of model fits
barPoints=[];
whichBars = [1:length(modelNamesAll)];

% modelNamesSort = modelNamesAll(whichBars);
for n=1:length(whichBars)
    for nr = numROIs
        tmp=barDataMeans(:, nr, :, :,n);
        tmp=(tmp(:));
        barPoints(:,n,nr)=tmp;
        tmp=tmp(~isnan(tmp));
        barMeans(n,nr)=mean(tmp);
        barStd(n,nr)=std(tmp);
        barSerr(n,nr)=std(tmp)/sqrt(length(tmp));
        CI95(n,:,nr) = tinv([0.025 0.975], length(tmp)-1);
        CI95(n,:,nr) =bsxfun(@times, barSerr(n,nr), CI95(n,:,nr));
    end
end
% figure;
sortBars = whichBars;
for nr = numROIs
%     subplot(length(numROIs),2,1+2*(nr-1));
    figure;
    subplot(1,2,1)
    bar(barMeans(sortBars,nr));
    hold on; errorbar(1:size(barMeans(sortBars,nr),1), barMeans(sortBars,nr), CI95(sortBars,1,nr), CI95(sortBars,2,nr), 'k.')
    axis square;
    xlim([0.5 length(modelNamesAll)+0.5])
    ylim([0 1])
    ax=gca;
    ax.XTickLabelRotation=90;
    ax.XTick = [1:length(modelNamesAll)];
    set(ax,'xticklabel',modelNamesAll(sortBars));
    ylabel('Median variance explained')
    title([mapNames{nr}])
    
    %Compute t-statistics and corresponding probabilities of pairwise
    %differences between model fits\
    for x=1:length(modelNamesAll)
        for y=1:length(modelNamesAll)
            [~,p,ci,stats] = ttest(barPoints(:,sortBars(x),nr), barPoints(:,sortBars(y),nr), 'tail', 'both');
            pvals(x,y,nr)=p;
            tMatrix(x,y,nr)=stats.tstat;
        end
    end
    tMatrix(tMatrix>200)=200;
    tMatrix(tMatrix<-200)=-200;
    
    %This matrix uses 10*10 pixel cells, because some viewing software
    %interpolates pixel edges
    tMatrixImg=zeros(size(tMatrix,1)*10, size(tMatrix,2)*10);
    for x=1:size(tMatrix,1)
        for y=1:size(tMatrix,1)
            tMatrixImg((x-1)*10+1:x*10, (y-1)*10+1:y*10)=tMatrix(x,y,nr);
        end
    end
    tMatrixImg(isnan(tMatrixImg))=0;
    %Add image of resulting t-statistics
    %figure;
%     subplot(length(numROIs),2,2*(nr));
    subplot(1,2,2)
    
    imagesc(tMatrixImg);
    clim = [-10 10];
    ax=gca;
    ax.XTickLabelRotation=90;
    ax.XTick = [10:10:length(modelNamesAll)*10];
    % set(ax,'xticklabel',modelNamesSort(sortBars))
    ax.YTick = [10:10:30];
    % set(ax,'yticklabel',modelNamesSort(sortBars))
    colormap(coolhotCmap([0], [128]));
    colorbar;
    axis image
    
    
end
% saveas(gcf, ['ModelComparisonWithinROI'], 'epsc');

%% Model Comparison all data Combined
thr= 0.2;

%
% modelNamesOriginal = modelNamesAll;
% modelNamesAll = modelNamesAll([3,4,1,2]);
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix
%Determine where to find data in structure
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1:length(modelNamesAll)
                % for variance explained
                targetDataOdd{whichSub, whichMap, whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Odd.ves'));
                targetDataEven{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Even.ves'));
            end
        end
    end
end

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1:length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataOdd{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataEven{whichSub, whichMap, whichHemi, modelN});
                    end
                    for oddEven=1:2
                        indices=max(barData{whichSub, whichMap, whichHemi, oddEven},[],2)>=0.2;
                        barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=median(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
                    end
                end
            end
        end
    end
end

%Compute and plot means and confidence intervals of model fits
barPoints=[];
whichBars = [1:length(modelNamesAll)];

% modelNamesSort = modelNamesAll(whichBars);
for n=1:length(whichBars)
    for nr = numROIs
        tmp=barDataMeans(:, :, :, :,n);
        tmp=(tmp(:));
        barPoints(:,n)=tmp;
        tmp=tmp(~isnan(tmp));
        barMeans(n)=mean(tmp);
        barStd(n)=std(tmp);
        barSerr(n)=std(tmp)/sqrt(length(tmp));
        CI95(n,:) = tinv([0.025 0.975], length(tmp)-1);
        CI95(n,:) =bsxfun(@times, barSerr(n), CI95(n,:));
    end
end
figure;
sortBars = whichBars;
subplot(1,2,1);
bar(barMeans(sortBars));
% hold on; errorbar(barMeans(sortBars), barMeans(sortBars), CI95(sortBars,1), CI95(sortBars,2), 'k.')
hold on; errorbar(1:size(barMeans,2), barMeans, CI95(:,1), CI95(:,2), 'k.')

axis square;
xlim([0.5 length(modelNamesAll)+0.5])
ylim([0 0.5])
ax=gca;
ax.XTickLabelRotation=90;
ax.XTick = [1:length(modelNamesAll)];
set(ax,'xticklabel',modelNamesAll(sortBars));
ylabel('Median variance explained')

%Compute t-statistics and corresponding probabilities of pairwise
%differences between model fits\
for x=1:length(modelNamesAll)
    for y=1:length(modelNamesAll)
        [~,p,ci,stats] = ttest(barPoints(:,sortBars(x)), barPoints(:,sortBars(y)), 'tail', 'both');
        pvals(x,y)=p;
        tMatrix(x,y)=stats.tstat;
    end
end
tMatrix(tMatrix>200)=200;
tMatrix(tMatrix<-200)=-200;

%This matrix uses 10*10 pixel cells, because some viewing software
%interpolates pixel edges
tMatrixImg=zeros(size(tMatrix,1)*10, size(tMatrix,2)*10);
for x=1:size(tMatrix,1)
    for y=1:size(tMatrix,1)
        tMatrixImg((x-1)*10+1:x*10, (y-1)*10+1:y*10)=tMatrix(x,y);
    end
end
tMatrixImg(isnan(tMatrixImg))=0;
%Add image of resulting t-statistics
%figure;
subplot(1,2,2);

imagesc(tMatrixImg);
ax=gca;
ax.XTickLabelRotation=90;
ax.XTick = [10:10:length(modelNamesAll)*10];
% set(ax,'xticklabel',modelNamesSort(sortBars))
ax.YTick = [10:10:30];
% set(ax,'yticklabel',modelNamesSort(sortBars))
colormap(coolhotCmap([0], [128]));
colorbar;
axis image
% end
saveas(gcf, ['ModelComparisonAllROIsCombined'], 'epsc');


% saveas(gcf, ['ModelComparisonWithinROI'], 'epsc');

%% Quick comparison between Odd/Even and Difficulty/Orientation
thr= 0.1;

%
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                targetDataOdd{whichSub, whichMap, whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Odd.x0s'));
                targetDataEven{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Even.x0s'));
                targetDataDifficulty{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Difficulty.x0s'));
                targetDataOrientation{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Orientation.x0s'));
                % Use original fits for index construction
                targetDataOddVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Odd.vesXval'));
                targetDataEvenVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Even.vesXval'));
                targetDataDifficultyVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Difficulty.vesXval'));
                targetDataOrientationVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Orientation.vesXval'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

for nm = 1:length(modelNamesAll)
    oePoints{nm}.points = [];
    doPoints{nm}.points = [];
end

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(subjectOrder)
    for whichMap=numROIs%length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataOdd{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataEven{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataDifficulty{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataOrientation{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 5}(:,modelN)=eval(targetDataOddVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 6}(:,modelN)=eval(targetDataEvenVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 7}(:,modelN)=eval(targetDataDifficultyVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 8}(:,modelN)=eval(targetDataOrientationVes{whichSub, whichMap, whichHemi, modelN});
                        %                     end
                        %                     for oddEven=1:2
                        %                         if oddEven == 1
                        %                             indices=max(barData{whichSub, whichMap, whichHemi, 5},[],2)>=thr;
                        indices=(barData{whichSub, whichMap, whichHemi, 5}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 6}(:,modelN)>=thr) & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=6  & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)>=1
                        oefullinds(modelN,whichSub,whichMap,whichHemi) = sum(indices);
%                         fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        %                         end
                        %                         barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=median(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
%                         oePoints{modelN}.points=[oePoints{modelN}.points;[barData{whichSub, whichMap, whichHemi, 1}(indices,modelN),barData{whichSub, whichMap, whichHemi, 2}(indices,modelN)]];
                        oePoints{modelN,whichSub,whichMap,whichHemi}.points=[barData{whichSub, whichMap, whichHemi, 1}(indices,modelN),barData{whichSub, whichMap, whichHemi, 2}(indices,modelN)];
                        % Make Orientation/Difficulty indices based on
                        % varExp and x-values below 10
                        indices=(barData{whichSub, whichMap, whichHemi, 7}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 8}(:,modelN)>=thr) & barData{whichSub, whichMap, whichHemi, 3}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=1;
                        dofullinds(modelN,whichSub,whichMap,whichHemi) = sum(indices);
%                         fullinds(2,modelN) = fullinds(2,modelN) + sum(indices);
                        doPoints{modelN,whichSub,whichMap,whichHemi}.points=[barData{whichSub, whichMap, whichHemi, 3}(indices,modelN),barData{whichSub, whichMap, whichHemi, 4}(indices,modelN)];
%                         doPoints{modelN}.points=[doPoints{modelN}.points;[barData{whichSub, whichMap, whichHemi, 3}(indices,modelN),barData{whichSub, whichMap, whichHemi, 4}(indices,modelN)]];
                    end
                end
            end
        end
    end
end

% uline = [1:11;1:11];

% figure;
for numModel = 1%length(modelNamesAll)
    for whichSub = 1:length(subjectOrder)
        for whichMap = 1:length(mapNames)
            for whichHemi = 1:2

    oddData = oePoints{numModel,whichSub,whichMap,whichHemi}.points(:,1);
    evenData = oePoints{numModel,whichSub,whichMap,whichHemi}.points(:,2);
    difficultyData = doPoints{numModel,whichSub,whichMap,whichHemi}.points(:,1);
    orientationData = doPoints{numModel,whichSub,whichMap,whichHemi}.points(:,2);
%     subplot(2,length(modelNamesAll),numModel)
%     subplot(2,2,numModel)
%     plot(oddData,evenData, '.')
%     hold on
%     plot(uline(1,:),uline(2,:))
%     axis([1 11 1 11])
%     axis square
%     title(['Odd vs Even, ', modelNamesAll{numModel}])
%     subplot(2,length(modelNamesAll),numModel+length(modelNamesAll))
%     subplot(2,2,numModel+2)
%     plot(difficultyData,orientationData, '.')
%     hold on
%     plot(uline(1,:),uline(2,:))
%     axis([1 11 1 11])
%     axis square
%     title(['Difficulty vs Orientation, ', modelNamesAll{numModel}])
    % numVoxels in model comparison corrected for upsampling to 1mm
    noe = oefullinds(numModel,whichSub,whichMap,whichHemi)./1.6^2;
    ndo = dofullinds(numModel,whichSub,whichMap,whichHemi)./1.6^2;
    
    tmp = corr(oddData,evenData, 'Type', 'Spearman');
    oer(numModel,whichSub,whichMap,whichHemi) = tmp;
    tmp = corr(difficultyData,orientationData, 'Type', 'Spearman');
    dor(numModel,whichSub,whichMap,whichHemi) = tmp;
    % correct for number of voxels and convert to p-val
%     oeip(numModel) = r2p(oer(numModel),noe);
%     doip(numModel) = r2p(dor(numModel),ndo);
    oeip(numModel,whichSub,whichMap,whichHemi) = r2p(oer(numModel,whichSub,whichMap,whichHemi),noe);
    doip(numModel,whichSub,whichMap,whichHemi) = r2p(dor(numModel,whichSub,whichMap,whichHemi),ndo);
               end
        end
    end 
end

% [oer,dor;oeip,doip]
[~,~,~,cor_oepvals] = fdr_bh(oeip, 0.05);
[~,~,~,cor_dopvals] = fdr_bh(doip, 0.05);

%% Plot stuff over distance
% showPlots(1) = all fits+upper+lower
% nsub = 1
% nps = 0;
% subjectOrder = {'dataAllSub'};
subjectOrder = {'DataS22'};
% subjectOrder = {'DataS01','DataS13','DataS14','DataS16', 'DataS19', 'DataS22'};

for nsub = 1:length(subjectOrder)
    % Left Hemi
    subName = subjectOrder{nsub}(end-2:end);
    clear data
    
    dataTmp = eval(subjectOrder{nsub});
    for nroi = 1:11;%length(mapNames)
        numMod = 1;
        for nm = 1
            data(nroi,numMod) = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).Progressive;
%             data(nroi+length(mapNames),numMod) = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).Progressive;
            numMod = numMod+1;
        end
    end
    veThresh = 0.20;
    showPlots = [1 0 0 0 0 0]; % Update to include the Odd/EvenćDiff/cOrie comparison plots
    binSteps = 21;
    meanThresh = [];
    figure;
    Hemi = {'Left'};
    data = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh, mapNames, subName, Hemi);
%     eval([subjectOrder{nsub}, 'edit = data']);
    for np = 1:length(mapNames)
        if np<=6
            if data(np).p<0.05
            npsips = npsips+1;
            else
            end
        end
        if np == 7 || np == 8
            if data(np).p<0.05
            npst = npst+1;
            else
            end
        end    
        if np >=9
            if data(np).p<0.05
            npsf = npsf+1;
            else
            end
        end
    end
    
    % Same for Right Hemi
    subName = subjectOrder{nsub}(end-2:end);
    clear data
    
    dataTmp = eval(subjectOrder{nsub});
    for nroi = 1:11;%length(mapNames)
        numMod = 1;
        for nm = 1
            data(nroi,numMod) = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).Progressive;
%             data(nroi+length(mapNames),numMod) = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).Progressive;
            numMod = numMod+1;
        end
    end
    veThresh = 0.20;
    showPlots = [1 0 0 0 0 0]; % Update to include the Odd/EvenćDiff/cOrie comparison plots
    binSteps = 21;
    meanThresh = [];
    figure;
    Hemi = {'Right'};
    data = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh, mapNames, subName, Hemi);
%     eval([subjectOrder{nsub}, 'edit = data']);
    for np = 1:length(mapNames)
        if np<=6
            if data(np).p<0.05
            npsips = npsips+1;
            else
            end
        end
        if np == 7 || np == 8
            if data(np).p<0.05
            npst = npst+1;
            else
            end
        end
        if np >= 9
            if data(np).p<0.05
            npsf = npsf+1;
            else
            end
        end
    end
end

% saveas(gcf, ['Load-versus-Distance'], 'epsc')


%% Plot Tuning Width stuff

% nsub = 1
subjectOrder = {'dataAllSub'}
% subjectOrder = {'DataS01','DataS13','DataS14','DataS16', 'DataS19', 'DataS22'};


for nsub = 1%:length(subjectOrder)
%     subName = subjectOrder{nsub}(end-2:end);
%     subName = 'DataAllSub';
    clear data
    
    dataTmp = eval(subjectOrder{1});
    for nroi = 1:11%length(mapNames)
        numMod = 1;
        for nm = 1
            data(nroi,numMod) = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).Progressive;
%             data(nroi+length(mapNames),numMod) = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).Progressive;
            numMod = numMod+1;
        end
    end
    veThresh = 0.1;
    showPlots = [0 0 0 1 1 0]; % Update to include the Odd/EvenćDiff/cOrie comparison plots
    binSteps = 21;
    meanThresh = [];
    figure;
    data = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh, mapNames, subName);
    eval([subjectOrder{nsub}, 'edit = data']);
end

%% Plot progressive fit + odd/even/orientation/difficulty lines
% showPlots(1) = all fits+upper+lower
% nsub = 1;
dts = {'Progressive', 'Orientation', 'Difficulty', 'Odd', 'Even'}%, 'RandomOrder'};
rois = [1:length(mapNames)];
% rois = 1
for nsub = [6]%:6]
%     if nsub == 4 || nsub == 6
%         dts = dts(1:end-1);
%     else
%         dts = {'Progressive', 'Orientation', 'Difficulty', 'Odd', 'Even'}%, 'RandomOrder'};
%     end
    subName = subjectOrder{nsub}(end-2:end);
    clear data
    dataTmp = eval(subjectOrder{nsub});
    for nroi = rois
        for ndts = 1:length(dts)
% Left            
%             data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).(dts{ndts}).meanDist;
%             data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).(dts{ndts}).ratio;
%             data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).(dts{ndts}).x0s;
%             data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).(dts{ndts}).ves;
%             data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).Left.(modelNamesAll{nm}).(dts{ndts}).sigmas;

% % Right
            data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).meanDist;
            data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).ratio;
            data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).x0s;
            data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).ves;
            data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).sigmas;
            
%             data(nroi+length(rois),ndts).sigmas = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).sigmas;


%             data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).meanDist;
%             data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).ratio;
%             data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).x0s;
%             data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).ves;
%             data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).sigmas;

        end
    end
    veThresh = 0.2;
    showPlots = [0 1 0 0 0 0]; % Update to include the Odd/EvenćDiff/cOrie comparison plots
    binSteps = 21; % 21 for 0.25 steps * range of 5 (6-1) + 1 for first step
    meanThresh = [];
%     figure;
    Hemi = {'Right'};
    data = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh, mapNames(rois), subName, Hemi);
    eval([subjectOrder{nsub}, 'edit = data']);
%     saveas(gcf, ['AllModelsOverDistance_', subName], 'epsc');
end

subjectOrder{nsub}
%%
for n = 1:11
    progressive(n) = data(n,1).p;
    orientation(n) = data(n,2).p;
    performance(n) = data(n,3).p;
    odd(n) = data(n,4).p;
    even(n) = data(n,5).p;
end

for n = 1:11
    allprog = (progressive<0.05) + (progressive<0.01) + (progressive<0.001) + (progressive<0.0001);
    allorient = (orientation<0.05) + (orientation<0.01) + (orientation<0.001) + (orientation<0.0001);
    allperf = (performance<0.05) + (performance<0.01) + (performance<0.001) + (performance<0.0001);
    allodd = (odd<0.05) + (odd<0.01) + (odd<0.001) + (odd<0.0001);
    alleven = (even<0.05) + (even<0.01) + (even<0.001) + (even<0.0001);
end

%% ANOVAs

%%
%ANOVAs
subOrder= subjectOrder;
hemispheres={'Left', 'Right'};
% backupNames = mapNames;
% backupModels = modelNamesAll;
modelNamesAll = modelNamesAll(1);
% tmp = mapNames([1,7,8,2,3,6,5,4,10,9,11])
% mapNames = tmp;
% mapNames = mapNames([7,8,1,2,[6:-1:3],[11:-1:9]])
% mapNames = mapNames([3,1:2,5:8,4,9:11])

x = 1;

VEs=nan([length(subjectOrder) length(mapNames) length(hemispheres)]);
SigmaMajor=VEs;
Q1D=VEs;
Q2D=VEs;
Q3D=VEs;
IQRD=VEs;
nVoxels=VEs;

% figure; hold on;
for n=1:length(subjectOrder)
    for m=1:length(mapNames)
        for whichHemi=1:length(hemispheres)
            if eval(char(strcat('isfield(', subjectOrder{n}, ', ''', mapNames{m}, ''') && isfield(', subjectOrder{n}, '.', mapNames{m}, ',''',hemispheres{whichHemi}, ''')' )))
                eval(char(strcat('data=', subjectOrder{n},'.', mapNames{m},'.', hemispheres{whichHemi},'.', modelNamesAll{1}  ,'.Progressive',';')));
                veIndices = data.ves > 0.1;
                if any(veIndices)
                    nVoxels(n,m,whichHemi)=sum(veIndices);
                end
                VEs(n,m,whichHemi)=mean(data.ves(veIndices));
                
                SigmaMajor(n,m,whichHemi)=mean(data.sigmas(veIndices));
                Q1D(n,m,whichHemi)=prctile(data.x0s(veIndices), 25);
                Q2D(n,m,whichHemi)=mean(data.x0s(veIndices));
                Q3D(n,m,whichHemi)=prctile(data.x0s(veIndices), 75);
                IQRD(n,m,whichHemi)=prctile(data.x0s(veIndices), 75)-prctile(data.x0s(veIndices), 25);
            end
        end
    end
end

%Setting up ANOVA structures
hemisphereGroups=cat(3, ones(size(VEs(:,:,1))), ones(size(VEs(:,:,1)))*2);
tmp=1:length(subjectOrder);
subjectLabels=cat(3, repmat(tmp(:), [1,length(mapNames)]), repmat(tmp(:), [1,length(mapNames)]));
mapLabels=cat(3, repmat(1:length(mapNames), [length(subjectOrder),1]), repmat(1:length(mapNames), [length(subjectOrder),1]));
mapLabels=mapLabels(:);
subjectLabels=subjectLabels(:);
hemisphereLabels=hemisphereGroups(:);
subjectLabels=subjectOrder(subjectLabels);
mapLabels=mapNames(mapLabels);
hemiTmp={'L', 'R'};
% hemisphereLabels=hemiTmp(hemisphereLabels);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME instead of ANOVA, apply to same list as all ANOVA's below

%create the data table, here I changed eccentricity from 0/1 to 1/2 here but I don't think it makes a difference
% fig2_data_lme=table(veANOVA,nominal(subjectANOVA),nominal(roiANOVA),nominal(eccANOVA+1),'VariableNames',{'ve', 'subject' ,'map', 'eccentricity'});

numSubs = repmat([1:length(subjectOrder)],[1,2*length(mapNames)]);
numMaps = repmat([1:length(mapNames)],[length(subjectOrder),2]);
% fig2_data_lme=table(VEs(:),nominal(numSubs(:)),nominal(numMaps(:)), nominal(hemisphereLabels(:)), nominal(Q2D(:)),'VariableNames',{'ve', 'subject' ,'map', 'hemisphere', 'exponent'});
% fig2_data_lme=table(VEs(:),nominal(numSubs(:)), nominal(numMaps(:)), nominal(Q2D(:)),'VariableNames',{'ve', 'subject', 'map', 'exponent'});
% fig2_data_lme=table(Exps(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'exponent', 'subject' ,'map'});

%pass the table into fitlme
%random factor is specified with (1|subject), 1 denotes intercept
%there are no interaction terms in this model just main effects, estimator is restricted maximum likelihood
% fig2_anova_lme= fitlme(fig2_data_lme,'ve~map+exponent+hemisphere+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
% fig2_anova_lme= fitlme(fig2_data_lme,'ve~1+map+(1|subject)','DummyVarCoding','effects','FitMethod','reml');
% fig2_anova_lme= fitlme(fig2_data_lme,'ve~map+exponent+(1|subject)','DummyVarCoding','reference','FitMethod','reml')
% fig2_anova_lme= fitlme(fig2_data_lme,'exponent~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')

%pass the lme results into anova to get df corrections for reporting combined main effects
% anova(fig2_anova_lme,'dfmethod','satterthwaite')

% nVox
fig2_data_lme=table(nVoxels(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'nVoxels', 'subject' ,'map'});
fig2_anova_lme= fitlme(fig2_data_lme,'nVoxels~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
anova(fig2_anova_lme,'dfmethod','satterthwaite')
% VE
fig2_data_lme=table(VEs(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'ve', 'subject' ,'map'});
fig2_anova_lme= fitlme(fig2_data_lme,'ve~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
anova(fig2_anova_lme,'dfmethod','satterthwaite')
% MeanD
fig2_data_lme=table(Q2D(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'meanDuration', 'subject' ,'map'});
fig2_anova_lme= fitlme(fig2_data_lme,'meanDuration~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
anova(fig2_anova_lme,'dfmethod','satterthwaite')
% IQRD
fig2_data_lme=table(IQRD(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'IQRD', 'subject' ,'map'});
fig2_anova_lme= fitlme(fig2_data_lme,'IQRD~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
anova(fig2_anova_lme,'dfmethod','satterthwaite')
% Major
fig2_data_lme=table(SigmaMajor(:),nominal(numSubs(:)),nominal(numMaps(:)),'VariableNames',{'SigmaMajor', 'subject' ,'map'});
fig2_anova_lme= fitlme(fig2_data_lme,'SigmaMajor~map+(1|subject)','DummyVarCoding','effects','FitMethod','reml')
anova(fig2_anova_lme,'dfmethod','satterthwaite')

%% ANOVA plots
%% Variance Explained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNIFICANT BETWEEN MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of VEs
% Hemi NOT significant
figure;
% [p, tmp, statsOut] = anovan(VEs(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
[p, tmp, statsOut] = anovan(VEs(:),{mapLabels subjectLabels }, 'varnames', {'map', 'subject'});
%Multiple comparison tests on ANOVA output
[results] =multcompare(statsOut, 'Dimension', [1]);

handle = gcf;
axObjs = handle.Children;
dataObjs = axObjs.Children;
for n=1:2:2*length(mapNames)-1; means((n+1)/2)=dataObjs(n+1).XData; CIs(((n+1)/2),:)=dataObjs(n).XData; end
means = fliplr(means);
CIs = flipud(CIs);
figure; plot(means, 'ok', 'MarkerFaceColor', 'k')
hold on;
for n = 1:length(mapNames)
    h = line([n n], [CIs(n,1) CIs(n,2)]);
    get(h)
    h.Color = [0 0 0];
    h.LineWidth = 2;
    set(h)
end
title(['Variance Explained'])
axis([0.5 length(mapNames)+0.5 0 1])
axis square
xticks([1:11])
xticklabels(mapNames)
xticklabel_rotate;
saveas(gcf, ['Variance Explained'], 'epsc');

%% Size in nVoxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNIFICANT BETWEEN MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of ROI voxel count
figure;
% Hemi NOT significant
% [p, tmp, statsOut] = anovan(nVoxels(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
[p, tmp, statsOut] = anovan(nVoxels(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});%Multiple comparison tests on ANOVA output
[results] =multcompare(statsOut, 'Dimension', [2]);
handle = gcf;
axObjs = handle.Children;
dataObjs = axObjs.Children;
for n=1:2:2*length(mapNames)-1; means((n+1)/2)=dataObjs(n+1).XData; CIs(((n+1)/2),:)=dataObjs(n).XData; end
means = fliplr(means);
CIs = flipud(CIs);
figure; plot(means, 'ok', 'MarkerFaceColor', 'k')
hold on;

for n = 1:length(mapNames)
    h = line([n n], [CIs(n,1) CIs(n,2)]);
    get(h)
    h.Color = [0 0 0];
    h.LineWidth = 2;
    set(h)
end
title(['number of Voxels'])
axis([0.5 length(mapNames)+0.5 0 1500])
axis square
xticks([1:11])
xticklabels(mapNames)
xticklabel_rotate;
saveas(gcf, ['number of Voxels'], 'epsc');

%% Sigma Majors
% Hemi NOT significant
% [p, tmp, statsOut] = anovan(SigmaMajor(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
figure;
[p, tmp, statsOut] = anovan(SigmaMajor(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
[results] =multcompare(statsOut, 'Dimension', [2]);
handle = gcf;
axObjs = handle.Children;
dataObjs = axObjs.Children;
for n=1:2:2*length(mapNames)-1; means((n+1)/2)=dataObjs(n+1).XData; CIs(((n+1)/2),:)=dataObjs(n).XData; end
means = fliplr(means);
CIs = flipud(CIs);
figure; plot(means, 'ok', 'MarkerFaceColor', 'k')
hold on;

for n = 1:length(mapNames)
    h = line([n n], [CIs(n,1) CIs(n,2)]);
    get(h)
    h.Color = [0 0 0];
    h.LineWidth = 2;
    set(h)
end
title(['Sigmas per Voxel'])
axis([0.5 length(mapNames)+0.5 0 30])
axis square
xticks([1:11])
xticklabels(mapNames)
xticklabel_rotate;
saveas(gcf, ['Sigma'], 'epsc');

%% Number Means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNIFICANT BETWEEN MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of duration mean
figure;
% Hemi NOT significant
% [p, tmp, statsOut] = anovan(Q2D(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
[p, tmp, statsOut] = anovan(Q2D(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
[results] =multcompare(statsOut, 'Dimension', [2]);
handle = gcf;
axObjs = handle.Children;
dataObjs = axObjs.Children;
for n=1:2:2*length(mapNames)-1; means((n+1)/2)=dataObjs(n+1).XData; CIs(((n+1)/2),:)=dataObjs(n).XData; end
means = fliplr(means);
CIs = flipud(CIs);
figure; plot(means, 'ok', 'MarkerFaceColor', 'k')
hold on;
for n = 1:length(mapNames)
    h = line([n n], [CIs(n,1) CIs(n,2)]);
    get(h)
    h.Color = [0 0 0];
    h.LineWidth = 2;
    set(h)
end
title(['Number means'])
axis([0.5 length(mapNames)+0.5 1 7])
axis square
xticks([1:11])
xticklabels(mapNames)
xticklabel_rotate;
saveas(gcf, ['Number Means'], 'epsc');
%Multiple comparison tests on ANOVA output

%% Number IQR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNIFICANT BETWEEN MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemi NOT significant
figure;
% [p, tmp, statsOut] = anovan(IQRD(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
[p, tmp, statsOut] = anovan(IQRD(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
%Multiple comparison tests on ANOVA output
[results] =multcompare(statsOut, 'Dimension', [2]);
handle = gcf;
axObjs = handle.Children;
dataObjs = axObjs.Children;
for n=1:2:2*length(mapNames)-1; means((n+1)/2)=dataObjs(n+1).XData; CIs(((n+1)/2),:)=dataObjs(n).XData; end
means = fliplr(means);
CIs = flipud(CIs);

figure; plot(means, 'ok', 'MarkerFaceColor', 'k')
hold on;
for n = 1:length(mapNames)
    h = line([n n], [CIs(n,1) CIs(n,2)]);
    get(h)
    h.Color = [0 0 0];
    h.LineWidth = 2;
    set(h)
end
title(['Number IQR'])
axis([0.5 length(mapNames)+0.5 0 5])
axis square
xticks([1:11])
xticklabels(mapNames)
xticklabel_rotate;
saveas(gcf, ['Number IQR'], 'epsc');


%% Quick comparison between Odd/Even and Difficulty/Orientation
thr= 0.1;

%
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2,3,5,6]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                targetDataNumerosity{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.RandomOrder.x0s'));
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
                targetDataNumerosityVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.RandomOrder.ves'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

for nm = 1:length(modelNamesAll)
    oePoints{nm}.points = [];
    doPoints{nm}.points = [];
end

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN});
                        %                     end
                        %                     for oddEven=1:2
                        %                         if oddEven == 1
                        %                             indices=max(barData{whichSub, whichMap, whichHemi, 5},[],2)>=thr;
                        indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr | barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr);% & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=10 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=10;
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        %                         end
                        %                         barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=median(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
                        oePoints{modelN}.points=[oePoints{modelN}.points;[barData{whichSub, whichMap, whichHemi, 1}(indices,modelN),barData{whichSub, whichMap, whichHemi, 2}(indices,modelN)]];
                    end
                end
            end
        end
    end
end


figure;
for numModel = 1%:length(modelNamesAll)
    oddData = oePoints{numModel}.points(:,1);
    evenData = oePoints{numModel}.points(:,2);
    subplot(2,ceil(length(modelNamesAll)/2),numModel)
    plot(oddData,evenData, '.')
    if numModel < 3
        axis([1 10 1 10])
    else
        axis([0 3 0 3])
    end
    title([modelNamesAll{numModel}])
    % numVoxels in model comparison corrected for upsampling to 1mm
    noe = fullinds(1,numModel)./1.6^2;
    
    
    tmp = corr(oddData,evenData, 'Type', 'Spearman');
    oer(numModel) = tmp;
    % correct for number of voxels and convert to p-val
%     oeip(numModel) = r2p(oer(numModel),noe);
    oeip(numModel) = r2p(oer(numModel),88);
    
    
end
suptitle(['Progressive versus RandomOrder'])

[oer;oeip]


%%
% r per map/hemisphere/subject, then compare all r values for sign and
% values (wilcoxon signed rank, paired difference). p-values for r, using
% nVoxels + upsampling correction.

% Odd vs Even
% Diff vs Orien
% Progr vs Random

%% Quick comparison between Odd/Even and Difficulty/Orientation
thr= 0.1;

%
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix rvals pvals
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2,3,5,6]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                targetDataNumerosity{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.RandomOrder.x0s'));
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
                targetDataNumerosityVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.RandomOrder.ves'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
    figure;
    for whichMap=1:length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:2%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN});

%                         indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr | barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr);% & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=10 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=10;
                        indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr) & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)>=1;
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        
                        test1 = barData{whichSub, whichMap, whichHemi, 1}(indices,modelN);
                        test2 = barData{whichSub, whichMap, whichHemi, 2}(indices,modelN);
                        try
                            r = corr(test1,test2, 'Type', 'Spearman');
                        catch
                            r= NaN;
                        end
                        rvals(whichSub,whichMap,whichHemi,modelN) = r;
                        pvals(whichSub,whichMap,whichHemi,modelN) = r2p(r,length(indices)./1.6^2);
%                         figure
                        subplot(length(hemispheres)+2,length(mapNames),whichMap+11*(modelN-1)+22*(whichHemi-1))
                        plot(test1,test2,'.')
                        hold on
                        plot([1:11],[1:11])
                        axis([1,11,1,11])
                        axis square
%                         shg
                    end
                end
            end
        end
    end
    suptitle([randomOrderSubs(whichSub)])
end

[~,~,~,cor_pvals] = fdr_bh(pvals, 0.05);
tmp=cor_pvals<0.05;
nps=sum(tmp(:));
%% VFM/Num vs VSTM

%% Quick comparison between Log2Lin/VFM/Numerosity
thr= 0.1;

%
numROIs = [1:length(mapNames)];
modelNamesAll = {'Log2Lin','Numerosity','VFM'};
m = 1:(length(modelNamesAll)); % whichModel

clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix rvals pvals
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
clear test1 test2 inds
indicesln = [];
indicesle = [];
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:2
            for modelN=1:length(modelNamesAll)
                % for variance explained
                if modelN == 3
                    targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat('sqrt((',randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s).^2 + (',randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.y0s).^2)'));
                else
                    targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                end
                
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
%     figure;
    for whichMap=1:length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1:3%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                    end
                    for modelN=1:3
                        if modelN == 1
                            inds1{whichSub,whichMap, whichHemi} = barData{whichSub, whichMap, whichHemi, 2}(:,1)>= thr & barData{whichSub, whichMap, whichHemi, 1}(:,1)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,1)<=8 &  barData{whichSub, whichMap, whichHemi, 1}(:,2)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,2)<=8;
                            inds2{whichSub,whichMap, whichHemi} = barData{whichSub, whichMap, whichHemi, 2}(:,1)>= thr & barData{whichSub, whichMap, whichHemi, 1}(:,1)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,1)<=8 &  barData{whichSub, whichMap, whichHemi, 1}(:,3)>=0 & barData{whichSub, whichMap, whichHemi, 1}(:,3)<=15;
                            indicesln= [indicesln; inds1{whichSub,whichMap, whichHemi}];
                            indicesle= [indicesle; inds2{whichSub,whichMap, whichHemi}];
                        else
                        end
                        
%                         fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        
                        test1{modelN} = barData{whichSub, whichMap, whichHemi, 1}(:,modelN);
                        if ~exist('test2', 'var')
                            for n = 1:length(modelNamesAll)
                                test2{n} = [];
                            end
%                             test2 = test1;
                        else
                        end
                          test2{modelN} = [test2{modelN}; test1{modelN}]; 
%                         test2 = barData{whichSub, whichMap, whichHemi, 2}(indices,modelN);
%                         r = corr(test1,test2, 'Type', 'Spearman');
%                         rvals(whichSub,whichMap,whichHemi,modelN) = r;
%                         pvals(whichSub,whichMap,whichHemi,modelN) = r2p(r,length(indices)./1.6^2);
%                         figure
%                         subplot(length(hemispheres)+2,length(mapNames),whichMap+11*(modelN-1)+22*(whichHemi-1))
%                         plot(test1,test2,'.')
%                         hold on
%                         plot([1:11],[1:11])
%                         axis([1,11,1,11])
%                         axis square
%                         shg
                    end
                end
            end
        end
    end
%     suptitle([randomOrderSubs(whichSub)])
end
%%

indicesln = logical(indicesln);
indicesle = logical(indicesle);
figure; 
subplot(1,2,1)
scatter(test2{1}, test2{2}, 'filled')

    r1 = corr(test2{1}(indicesln), test2{2}(indicesln), 'Type', 'Spearman');
    p1 = r2p(r1,sum(indicesln)./1.6^2);
    p1 = r2p(r1,44);
xlabel('VSTM Load')
ylabel('Numerosity')
title([r])
axis square


subplot(1,2,2)
scatter(test2{1}, test2{3}, 'filled')
    r2 = corr(test2{1}(indicesle), test2{3}(indicesle), 'Type', 'Spearman');    
    p2 = r2p(r2,sum(indicesle)./1.6^2);
    p2 = r2p(r2,44);
    
xlabel('VSTM Load')
ylabel('Eccenctricity')
title([r])
axis square



%% Test how many maps are significantly increasing per participant.

% showPlots(1) = all fits+upper+lower
% nsub = 1;
% subjectOrder = {'DataS14'};
subjectOrder = {'DataS01','DataS13','DataS14','DataS16', 'DataS19', 'DataS22'};
dts = {'Progressive', 'Orientation', 'Difficulty', 'Odd', 'Even'}%, 'RandomOrder'};
rois = [1:length(mapNames)];
Hemi = {'Left', 'Right'}
% rois = 1
for nh = 1:2
for nsub = [1:6]
%     if nsub == 4 || nsub == 6
%         dts = dts(1:end-1);
%     else
%         dts = {'Progressive', 'Orientation', 'Difficulty', 'Odd', 'Even'}%, 'RandomOrder'};
%     end
    subName = subjectOrder{nsub}(end-2:end);
    clear data
    dataTmp = eval(subjectOrder{nsub});
    for nroi = rois
        for ndts = 1:length(dts)
            data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).(Hemi{nh}).(modelNamesAll{nm}).(dts{ndts}).meanDist;
            data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).(Hemi{nh}).(modelNamesAll{nm}).(dts{ndts}).ratio;
            data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).(Hemi{nh}).(modelNamesAll{nm}).(dts{ndts}).x0s;
            data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).(Hemi{nh}).(modelNamesAll{nm}).(dts{ndts}).ves;
            data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).(Hemi{nh}).(modelNamesAll{nm}).(dts{ndts}).sigmas;
%             data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).meanDist;
%             data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).ratio;
%             data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).x0s;
%             data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).ves;
%             data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).sigmas;
            
%             data(nroi+length(rois),ndts).sigmas = dataTmp.(mapNames{nroi}).Right.(modelNamesAll{nm}).(dts{ndts}).sigmas;


%             data(nroi,ndts).meanDist = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).meanDist;
%             data(nroi,ndts).ratio = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).ratio;
%             data(nroi,ndts).x0s = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).x0s;
%             data(nroi,ndts).ves = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).ves;
%             data(nroi,ndts).sigmas = dataTmp.(mapNames{nroi}).Both.(modelNamesAll{nm}).(dts{ndts}).sigmas;

        end
    end
    veThresh = 0.2;
    showPlots = [0 1 0 0 0 0]; % Update to include the Odd/EvenćDiff/cOrie comparison plots
    binSteps = 21; % 21 for 0.25 steps * range of 5 (6-1) + 1 for first step
    meanThresh = [];
%     figure;
    data = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh, mapNames(rois), subName, Hemi(nh));
    eval([subjectOrder{nsub}, Hemi{nh}, ' = data']);
%     saveas(gcf, ['AllModelsOverDistance_', subName], 'epsc');
end
end

%% Load vs Numerosity
thr= 0.1;

%
modelNamesAll = {'Log2Lin', 'Numerosity'};
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix rvals pvals
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2,4]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:length(hemispheres)
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                targetDataNumerosity{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+1}, '.Progressive.x0s'));
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
                targetDataNumerosityVes{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+1}, '.Progressive.ves'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
    figure;
    for whichMap=1:length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:length(hemispheres)
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:2%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN+1});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN+1});

                        indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)>0 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=7 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr);
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        
                        test1 = barData{whichSub, whichMap, whichHemi, 1}(indices,modelN);
                        test2 = barData{whichSub, whichMap, whichHemi, 2}(indices,modelN);
                        try
                            r = corr(test1,test2, 'Type', 'Spearman');
                        catch
                            r = NaN;
                        end
                        indset(whichSub,whichMap,whichHemi,modelN) = sum(indices);
                        rvals(whichSub,whichMap,whichHemi,modelN) = r;
                        pvals(whichSub,whichMap,whichHemi,modelN) = r2p(r,sum(indices)./1.6^2);
%                         figure
                        subplot(length(hemispheres)+2,length(mapNames),whichMap+11*(modelN-1)+22*(whichHemi-1))
                        plot(test1,test2,'.')
                        hold on
                        plot([1:11],[1:11])
                        axis([1,11,1,11])
                        axis square
%                         shg
                    end
                end
            end
        end
    end
    suptitle([randomOrderSubs(whichSub)])
end

[~,~,~,cor_pvals] = fdr_bh(pvals, 0.05);

%% Load vs VFM
thr= 0.1;
% modelNamesAll = {'Log2Lin', 'Numerosity', 'VFM'};
modelNamesAll = {'Log2Lin', 'VFM'};
% modelNamesAll = {'Log2Lin', 'Numerosity'};

%
% hemispheres = {'Both'};
hemispheres = {'Left', 'Right'};
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix rvals pvals
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2,3,5,6]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:length(hemispheres)
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                targetDataNumerosity{whichSub, whichMap, whichHemi, modelN+1}=char(strcat('sqrt((',randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+1}, '.Progressive.x0s).^2 + (',randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+1}, '.Progressive.y0s).^2)'));
%                 targetDataNumerosity{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+2}, '.Progressive.x0s'));
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
                targetDataNumerosityVes{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+1}, '.Progressive.ves'));
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
    figure;
    for whichMap=1:length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:length(hemispheres)
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:2%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN+1});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN+1});
                        
%                         if whichMap == 9
%                             x = 1;
%                         end
                        
                        indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)>0 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=15 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr) 
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        
                        test1 = barData{whichSub, whichMap, whichHemi, 1}(indices,modelN);
                        test2 = barData{whichSub, whichMap, whichHemi, 2}(indices,modelN);
                        try
                            r = corr(test1,test2, 'Type', 'Spearman');
                        catch
                            r= NaN;
                        end
                        indset(whichSub,whichMap,whichHemi,modelN) = sum(indices);
                        rvals(whichSub,whichMap,whichHemi,modelN) = r;
                        pvals(whichSub,whichMap,whichHemi,modelN) = r2p(r,sum(indices)./1.6^2);
%                         figure
                        subplot(length(hemispheres)+2,length(mapNames),whichMap+11*(modelN-1)+22*(whichHemi-1))
                        plot(test1,test2,'.')
                        hold on
                        plot([1:11],[1:11])
                        axis([1,11,1,11])
                        axis square
%                         shg
                    end
                end
            end
        end
    end
    suptitle([randomOrderSubs(whichSub)])
end

[~,~,~,cor_pvals] = fdr_bh(pvals, 0.05);
%%
%% Load vs Reverse
thr= 0.1;

modelNamesAll = {'Log2Lin', 'Numerosity', 'VFM', 'Reverse'};
%
% hemispheres = {'Both'};
hemispheres = {'Left', 'Right'};
numROIs = [1:length(mapNames)];
m = 1:(length(modelNamesAll)); % whichModel
clear barDataMeans barPoints barMeans barStd barSerr CI95 barData targetDataOdd targetDataEven indices tMatrix rvals pvals
clear oer dor oeip doip fullinds oePoints doPoints oddData evenData orientationData difficultyData
%Determine where to find data in structure
randomOrderSubs = subjectOrder([1,2,6]);
for whichSub= 1:length(randomOrderSubs)
    for whichMap=numROIs%length(mapNames)
        for whichHemi=1:length(hemispheres)
            for modelN=1%:length(modelNamesAll)
                % for variance explained
                if whichSub == 1 || whichSub == 3
                    targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.Log2Lin.SmallDifficulty.x0s'));
                elseif whichSub == 2
                    targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.Log2Lin.Difficulty.x0s'));
                else
                end
                        
%                         targetDataProgressive{whichSub, whichMap, whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.x0s'));
                %                 targetDataNumerosity{whichSub, whichMap, whichHemi, modelN+1}=char(strcat('exp(', randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN+3}, '.Progressive.x0s)'));
%                 if whichSub == 3
                    targetDataNumerosity{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.Reverse.Log2Lin.x0s'));
%                 else
%                 end
                
                % Use original fits for index construction
                targetDataProgressiveVes{whichSub, whichMap,whichHemi, modelN}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.Progressive.ves'));
                if whichSub == 3
                    targetDataNumerosityVes{whichSub, whichMap,whichHemi, modelN+1}=char(strcat(randomOrderSubs{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.Reverse.Log2Lin.ves'));
                else
                end
            end
        end
    end
end

fullinds = zeros(2,length(modelNamesAll));

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(randomOrderSubs), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(randomOrderSubs)
    figure;
    for whichMap=1:length(mapNames)
        if isfield(eval(randomOrderSubs{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:length(hemispheres)
                if isfield(eval([char(randomOrderSubs{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:2%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN+1});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN+1});
                        
%                         if whichMap == 9
%                             x = 1;
%                         end
                        
                        indices=(barData{whichSub, whichMap, whichHemi, 3}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)>=1 & barData{whichSub, whichMap, whichHemi, 1}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)>1 & barData{whichSub, whichMap, whichHemi, 2}(:,modelN)<=6 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN)>=thr);
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        
                        test1 = barData{whichSub, whichMap, whichHemi, 1}(indices,modelN);
                        test2 = barData{whichSub, whichMap, whichHemi, 2}(indices,modelN);
                        try
                            r = corr(test1,test2, 'Type', 'Spearman');
                        catch
                            r= NaN;
                        end
                        indset(whichSub,whichMap,whichHemi,modelN) = sum(indices);
                        rvals(whichSub,whichMap,whichHemi,modelN) = r;
                        pvals(whichSub,whichMap,whichHemi,modelN) = r2p(r,sum(indices)./1.6^2);
%                         figure
                        subplot(length(hemispheres)+2,length(mapNames),whichMap+11*(modelN-1)+22*(whichHemi-1))
                        plot(test1,test2,'.')
                        hold on
                        plot([1:11],[1:11])
                        axis([1,6,1,6])
                        axis square
%                         shg
                    end
                end
            end
        end
    end
    suptitle([randomOrderSubs(whichSub)])
end

[~,~,~,cor_pvals] = fdr_bh(pvals, 0.05);
tmp = cor_pvals<0.05;
nps = sum(tmp(:));
% modelNamesAll = {'Log2Lin', 'Numerosity'}
%%
progdat = [];
vfmdat = [];
for ns = 1:2
    for nm = 1:11
        progdat = [progdat; [barData{ns,nm,1,1}((barData{ns,nm,1,3}>thr & barData{ns,nm,1,4}>thr),1)]];
        vfmdat = [vfmdat; [barData{ns,nm,1,2}((barData{ns,nm,1,3}>thr & barData{ns,nm,1,4}>thr),1)]];
    end
end

numinds = length(progdat);

r = corr(progdat, vfmdat, 'Type', 'Spearman');
p = r2p(r,numinds./1.6^2)



