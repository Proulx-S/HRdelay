function [resPall,resAll] = permuteLabels(p,resAll)
repoPath = p.paths.repo.in;
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir = 'd';
subjList = p.meta.subjList;
sessList = {'sess1' 'sess2'};

% disp('Permuting Labels: Permuting')
if ~exist('resAll','var') || isempty(resAll)
    resAll = cell(length(subjList),1);
    loadFlag = 1;
else
    loadFlag = 0;
end
resPall = cell(length(subjList),1);
for subjInd = 1:length(subjList)
    subj = subjList{subjInd};
    %Load
    if loadFlag
        load(fullfile(funPath,inDir,[subj '.mat']),'res')
        resAll{subjInd} = res;
    else
        res = resAll{subjInd};
    end
    resP = res;
    for sessInd = 1:length(sessList)
        sess = sessList{sessInd};
        for rep = 1:size(res.(sess).sin,4)
            %Permute
            perm = randperm(3);
            resP.(sess).sin(:,:,:,rep,:) = res.(sess).sin(:,:,:,rep,perm);
            resP.(sess).sinBase(:,:,:,rep,:) = res.(sess).sinBase(:,:,:,rep,perm);
            resP.(sess).hr(:,:,:,rep,:,:) = res.(sess).hr(:,:,:,rep,perm,:);
            resP.(sess).hrBase(:,:,:,rep,:,:) = res.(sess).hrBase(:,:,:,rep,perm,:);
        end
    end
    resPall{subjInd} = resP;
    %Save
%     save(fullfile(funPath,outDir,[subj '.mat']),'resP','-append')
end
% disp('Permuting Labels: Done')