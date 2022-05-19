function [respPall,respAll] = permuteLabels(p,respAll)
% repoPath = p.paths.repo.in;
        funPath = fullfile(p.dataPath.V1,'resp');
%             inDir  = 'd';
%             outDir = 'd';
subjList = p.meta.subjList;
sessList = {'sess1' 'sess2'};

% disp('Permuting Labels: Permuting')
if ~exist('resAll','var') || isempty(respAll)
    respAll = cell(length(subjList),1);
    loadFlag = 1;
else
    loadFlag = 0;
end
respPall = cell(length(subjList),1);
for subjInd = 1:length(subjList)
    subj = subjList{subjInd};
    %Load
    if loadFlag
        load(fullfile(funPath,[subj '.mat']),'resp')
        respAll{subjInd} = resp;
    else
        resp = respAll{subjInd};
    end
    respP = resp;
    for sessInd = 1:length(sessList)
        sess = sessList{sessInd};
        for rep = 1:size(resp.(sess).sin,4)
            %Permute
            perm = randperm(3);
            respP.(sess).sin(:,:,:,rep,:) = resp.(sess).sin(:,:,:,rep,perm);
            respP.(sess).sinBase(:,:,:,rep,:) = resp.(sess).sinBase(:,:,:,rep,perm);
            respP.(sess).hr(:,:,:,rep,:,:) = resp.(sess).hr(:,:,:,rep,perm,:);
            respP.(sess).hrBase(:,:,:,rep,:,:) = resp.(sess).hrBase(:,:,:,rep,perm,:);
        end
    end
    respPall{subjInd} = respP;
    %Save
%     save(fullfile(funPath,outDir,[subj '.mat']),'resP','-append')
end
% disp('Permuting Labels: Done')