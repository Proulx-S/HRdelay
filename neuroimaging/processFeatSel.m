function processFeatSel(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,2)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
%     for sessInd = 1:2
%         sess = ['sess' num2str(sessInd)];
%         if isfield(res.(sess),'featSel')
%             res.(sess) = rmfield(res.(sess),'featSel');
%         end
%     end
    dAll{subjInd} = res;
end
d = dAll; clear dAll
sessList = fields(d{1});


%% Reorganize
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        sess = ['sess' num2str(sessInd)];
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
    end
end
d = dP; clear dP

%% Precompute flattened voxel ecc distribution on fov
if p.featSel.fov.doIt && strcmp(p.featSel.fov.threshMethod,'empirical')
    disp('Flattening: computing')
    d = flattenEccDist(d,p,1);
    disp('Flattening: done')
    p.featSel.fov.empirical.padFac = 1.2;
    d = prepareDelayFovContour(d,p);
    %                                           1    2    3    4    5    6
    p.featSel.fov.empirical.smList           = [0.01  0.15  0.15  0.15  0.15  0.15
                                                0.01  0.15  0.15  0.15  0.15  0.15]; % ecc
    p.featSel.fov.empirical.levelList        = [0.50  0.50  0.50  0.50  0.50  0.50
                                                0.50  0.50  0.50  0.50  0.50  0.50]; % 0:1
    p.featSel.fov.empirical.mergeRadiusList  = [0.60  0.60  0.60  0.60  0.60  0.60
                                                0.60  0.60  0.60  0.60  0.60  0.60]; % ecc
    p.featSel.fov.empirical.marginRadiusList = [0.40  0.40  0.40  0.40  0.40  0.40
                                                0.40  0.40  0.40  0.40  0.40  0.40]; % ecc
    p.featSel.fov.empirical.contIndList1     = {[1 4 5 11 12 16 20 45  2 8 7 3 37 38 9 39 26 29 42 6 14 24] [inf] [inf] [inf] [inf] [inf]
                                                [1 4 5 12 18 26 48  9 21 35 3 7 14 2 29 33 8 32 6] [inf] [inf] [inf] [inf] [inf]}; % contour indices, positive values->select voxels inside contour; negative values->select voxels outside contour
    p.featSel.fov.empirical.contIndList2     = {[inf] [inf] [inf] [inf] [inf] [inf]
                                                [inf] [inf] [inf] [inf] [inf] [inf]}; % contour indices, positive values->select voxels inside contour; negative values->select voxels outside contour

end

%% Feature selection
featSel = cell(size(d));
f = cell(size(d));
% disp('computing feature selection stats')
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        p.subjInd = subjInd;
        p.sessInd = sessInd;
        disp(['subj' num2str(subjInd) '; sess' num2str(sessInd)])
        [featSel{subjInd,sessInd},f{subjInd,sessInd}] = getFeatSel(d{subjInd,sessInd},p);
    end
end

fAll = cell(size(f{1}));
supTitleList = {'Contour Definition' 'Contour Processing' 'Final Contour' 'Contour Masking'};
for fAllInd = 1:size(f{1},2)
    if fAllInd~=4
        fAll{fAllInd} = figure('WindowStyle','docked');
        [ha, pos] = tight_subplot(size(d,2), size(d,1), 0, 0.1, 0); delete(ha);
        for subjInd = 1:size(d,1)
            for sessInd = 1:size(d,2)
                ax = copyobj(f{subjInd,sessInd}(fAllInd).Children,fAll{fAllInd});
                ax.Position = pos{(sessInd-1)*size(d,1)+subjInd};
                ax.Colormap = f{subjInd,sessInd}(fAllInd).Children.Colormap;
                delete(f{subjInd,sessInd}(fAllInd).Children);
            end
        end
        suptitle(supTitleList{fAllInd})
    end
end


indInX = cell(1,size(d,2));
dX = cell(1,size(d,2));
fieldList = fields(d{subjInd,sessInd});
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        if subjInd==1
            indInX{sessInd} = featSel{subjInd,sessInd}.featSeq.featIndIn;
            dX{sessInd} = d{subjInd,sessInd};
        else
            indInX{sessInd} = cat(1,indInX{sessInd},featSel{subjInd,sessInd}.featSeq.featIndIn);
            for fieldInd = 1:length(fieldList)
                if isnumeric(d{subjInd,sessInd}.(fieldList{fieldInd}))...
                        && ~strcmp(fieldList{fieldInd},'sinDesign')...
                        && ~strcmp(fieldList{fieldInd},'hrDesign')
                    dX{sessInd}.(fieldList{fieldInd}) = cat(1,mean(dX{sessInd}.(fieldList{fieldInd}),4),mean(d{subjInd,sessInd}.(fieldList{fieldInd}),4));
                elseif isstruct(d{subjInd,sessInd}.(fieldList{fieldInd}))...
                        && ~strcmp(fieldList{fieldInd},'featSel')
                    fieldList2 = fields(d{subjInd,sessInd}.(fieldList{fieldInd}));
                    for fieldInd2 = 1:length(fieldList2)
                        if isnumeric(d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}))...
                                || islogical(d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}))
                            dX{sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}) = ...
                                cat(1,dX{sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}),d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}));    
                        end
                    end
                end
            end
        end
    end
end

sessInd = 1;
plotVoxOnFoV(dX{sessInd},p,true(size(dX{sessInd}.sin,1),1))
ax = gca;
RLim = ax.RLim;
title('allVox')
for featInd = 1:length(featSel{subjInd,sessInd}.featSeq.featSelList)
    plotVoxOnFoV(dX{sessInd},p,indInX{sessInd}(:,featInd,1))
    title(featSel{subjInd,sessInd}.featSeq.featSelList{featInd})
    ax = gca;
    ax.RLim = RLim;
end

featSel{subjInd,sessInd}.featSeq.featSelList'
featInd = [3 4];
plotVoxOnFoV(dX{sessInd},p,all(indInX{sessInd}(:,featInd,1),2))
ax = gca;
R = [ax.Children(3).RData ax.Children(4).RData];
figure('WindowStyle','docked');
R = exp(R)-1;
hist(R,100)




featInd = [0];
condPairInd = 1;
if featInd==0
    ind = true(size(featSel{subjInd,sessInd}.featSeq.featIndIn,1),1);
else
    ind = all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,featInd,condPairInd),2);
end
p.featSel.fov.threshVal = [];
plotVoxOnFoV(dX{1},p,ind)

% fac = 1;
% fac = (2*pi);
% fac = 1/(2*pi);
% fac = (2*2*pi);
fac = 1/(2*2*pi);

R = dX{sessInd}.voxProp.ecc;
Rp = cdf(nonparamDistFit(R),R);
v = linspace(min(R),max(R)*fac,length(unique(R)));%R
x = linspace(1/length(unique(R)),1,length(unique(R)));%Rp
xq = Rp;%Rp
vq = interp1(x,v,xq);%R
R2 = vq;
Rp2 = xq;
figure('WindowStyle','docked');
[~,uInd,~] = unique(R);
plot(R(uInd),Rp(uInd),'.'); hold on
[~,uInd,~] = unique(R2);
plot(R2(uInd),Rp2(uInd),'.'); hold on
nonparamDistFit(R2,0)

d2 = dX{sessInd};
d2.voxProp.ecc = R2;
plotVoxOnFoV(d2,p,ind)





%% Save
disp('saving feature selection')
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
save(fullfilename,'featSel')
disp(['saved to: ' fullfilename])
