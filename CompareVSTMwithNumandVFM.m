mapNames = ["VMIPS0","VMIPS1","VMIPS2","VMIPS3","VMIPS4","VMIPS5","VMT1","VMT2","VMF1","VMF2","VMF3"];
mapNames = ["VMIPS0","VMT1","VMT2","VMIPS1","VMIPS2","VMIPS3","VMIPS4","VMIPS5","VMF3","VMF2","VMF1"];
% subjectOrder = {'DataS01','DataS13','DataS14','DataS16','DataS19'};
subjectOrder = {'DataS01','DataS13','DataS14','DataS16', 'DataS19', 'DataS22'};
hemispheres = {'Left', 'Right'};
% modelNamesAll =     {'Lin','Log2Lin', 'LinMonotonic', 'LogMonotonic'};
modelNamesAll = {'Log2Lin', 'VFM', 'Numerosity'};


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
            % for variance explained
            targetDataProgressive{whichSub, whichMap, whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{1}, '.Progressive.x0s'));
            targetDataVFMx{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{2}, '.Progressive.x0s'));
            targetDataVFMy{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{2}, '.Progressive.y0s'));
            targetDataNumerosity{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{3}, '.Progressive.x0s'));
            % Use original fits for index construction
            targetDataProgressiveVes{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{1}, '.Progressive.ves'));
            targetDataVFMVes{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{2}, '.Progressive.ves'));
            targetDataNumerosityVes{whichSub, whichMap,whichHemi, 1}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{3}, '.Progressive.ves'));
        end
    end
end


% clear rs
% clear ps

rs = [];
ps = [];
for mm = 1:11
fullinds = zeros(3,1);

for nm = 1:length(modelNamesAll)
    VFMPoints{nm}.points = [];
    NumPoints{nm}.points = [];
end


%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub= 1:length(subjectOrder)
    for whichMap=mm%numROIs%length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1%:2%length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataProgressive{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataVFMx{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 3}(:,modelN)=eval(targetDataVFMy{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 4}(:,modelN)=eval(targetDataNumerosity{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 5}(:,modelN)=eval(targetDataProgressiveVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 6}(:,modelN)=eval(targetDataVFMVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 7}(:,modelN)=eval(targetDataNumerosityVes{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 8}(:,modelN)=sqrt(barData{whichSub, whichMap, whichHemi,2}.^2+barData{whichSub, whichMap, whichHemi,3}.^2);
                        
                        
                        indices=barData{whichSub, whichMap, whichHemi, 5}(:,modelN)>=thr;
                        VFMindices=(barData{whichSub, whichMap, whichHemi, 6}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 5}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 1}(:,modelN) < 8 & barData{whichSub, whichMap, whichHemi, 1}(:,modelN) > 1 & barData{whichSub, whichMap, whichHemi, 8}(:,modelN) < 20);
                        Numindices=(barData{whichSub, whichMap, whichHemi, 7}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 5}(:,modelN)>=thr & barData{whichSub, whichMap, whichHemi, 1}(:,modelN) < 8 & barData{whichSub, whichMap, whichHemi, 4}(:,modelN) < 8);
                        fullinds(1,modelN) = fullinds(1,modelN) + sum(indices);
                        fullinds(2,modelN) = fullinds(2,modelN) + sum(VFMindices);
                        fullinds(3,modelN) = fullinds(3,modelN) + sum(Numindices);
                        
                        VFMPoints{modelN}.points=[VFMPoints{modelN}.points;[barData{whichSub, whichMap, whichHemi, 1}(VFMindices,modelN),barData{whichSub, whichMap, whichHemi, 8}(VFMindices,modelN)]];
                        NumPoints{modelN}.points=[NumPoints{modelN}.points;[barData{whichSub, whichMap, whichHemi, 1}(Numindices,modelN),barData{whichSub, whichMap, whichHemi, 4}(Numindices,modelN)]];
                        % Make Orientation/Difficulty indices based on
                        % varExp and x-values below 10
                    end
                end
            end
        end
    end
end

uline = [1:11;1:11];
numModel = 1
figure;
VSTMVFMData = VFMPoints{numModel}.points(:,1);
VFMData = VFMPoints{numModel}.points(:,2);
VSTMNumData = NumPoints{numModel}.points(:,1);
NumData = NumPoints{numModel}.points(:,2);
%     subplot(2,length(modelNamesAll),numModel)
subplot(1,2,1)
plot(VSTMVFMData,VFMData, '.')
hold on
plot(uline(1,:),uline(2,:))
% axis([1 11 1 11])
axis square
title(['VSTM vs Eccentricity'])
%     subplot(2,length(modelNamesAll),numModel+length(modelNamesAll))
subplot(1,2,2)
plot(VSTMNumData,NumData, '.')
hold on
plot(uline(1,:),uline(2,:))
% axis([1 11 1 11])
axis square
title(['VSTM vs Numerosity'])
% numVoxels in model comparison corrected for upsampling to 1mm
nvfm = fullinds(2,numModel)./1.6^2;
nnum = fullinds(3,numModel)./1.6^2;

tmp = corr(VSTMVFMData,VFMData, 'Type', 'Spearman');
vfmr = tmp;
tmp = corr(VSTMNumData,NumData, 'Type', 'Spearman');
numr = tmp;
% correct for number of voxels and convert to p-val
%     oeip(numModel) = r2p(oer(numModel),noe);
%     doip(numModel) = r2p(dor(numModel),ndo);
vfmp = r2p(vfmr,nvfm);
nump = r2p(numr,nnum);


rs = [rs;[vfmr,numr]]
ps = [ps;[vfmp,nump]]
end
