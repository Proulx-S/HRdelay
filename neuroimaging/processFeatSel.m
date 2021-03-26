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


%% Feature selection
featSel = cell(size(d));
disp('computing feature selection stats')
for sessInd = 1:numel(d)
    disp(['sess' num2str(sessInd) '/' num2str(numel(d))])
    [subjInd,sessInd] = ind2sub(size(d),sessInd);
    featSel{subjInd,sessInd} = getFeatSel(d{subjInd,sessInd},p);
end
disp('done')


warning('finish coding that or comment it out')
keyboard
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
tmpR = [ax.Children(3).RData ax.Children(4).RData];
figure('WindowStyle','docked');
tmpR = exp(tmpR)-1;
hist(tmpR,100)


%% Save
disp('saving feature selection')
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
save(fullfilename,'featSel')
disp(['saved to: ' fullfilename])
