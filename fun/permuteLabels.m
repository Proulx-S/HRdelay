function [respP,resp] = permuteLabels(p,resp)
% repoPath = p.paths.repo.in;
        funPath = fullfile(p.dataPath.V1,'resp');
%             inDir  = 'd';
%             outDir = 'd';
subjList = p.meta.subjList;
sessList = {'sess1' 'sess2'};

%% Load data if not provided as input respAll
if ~exist('respAll','var') || isempty(resp)
    respAll = cell(length(subjList),1);
    for subjInd = 1:length(subjList)
        subj = subjList{subjInd};
        %Load
        load(fullfile(funPath,[subj '.mat']),'resp')
        respAll{subjInd} = resp;
    end
end
[resp, info] = reorgData(p,respAll); clear respAll

%% Permute session by session, run repetition (triplet) by run repetition 
respP = resp;
for sessInd = 1:numel(resp)
    for repInd = 1:size(resp{sessInd}.sin,4)
        %Permute
        perm = randperm(3);
        respP{sessInd}.sin(:,:,:,repInd,:) = resp{sessInd}.sin(:,:,:,repInd,perm);
        respP{sessInd}.sinBase(:,:,:,repInd,:) = resp{sessInd}.sinBase(:,:,:,repInd,perm);
        respP{sessInd}.hr(:,:,:,repInd,:,:) = resp{sessInd}.hr(:,:,:,repInd,perm,:);
        respP{sessInd}.hrBase(:,:,:,repInd,:,:) = resp{sessInd}.hrBase(:,:,:,repInd,perm,:);
    end
end



% % disp('Permuting Labels: Permuting')
% if ~exist('respAll','var') || isempty(resp)
%     resp = cell(length(subjList),1);
%     loadFlag = 1;
% else
%     loadFlag = 0;
% end
% respP = cell(length(subjList),1);
% for subjInd = 1:length(subjList)
%     subj = subjList{subjInd};
%     %Load
%     if loadFlag
%         load(fullfile(funPath,[subj '.mat']),'resp')
%         resp{subjInd} = resp;
%     else
%         resp = resp{subjInd};
%     end
%     respP = resp;
%     for sessInd = 1:length(sessList)
%         sess = sessList{sessInd};
%         for rep = 1:size(resp.(sess).sin,4)
%             %Permute
%             perm = randperm(3);
%             respP.(sess).sin(:,:,:,rep,:) = resp.(sess).sin(:,:,:,rep,perm);
%             respP.(sess).sinBase(:,:,:,rep,:) = resp.(sess).sinBase(:,:,:,rep,perm);
%             respP.(sess).hr(:,:,:,rep,:,:) = resp.(sess).hr(:,:,:,rep,perm,:);
%             respP.(sess).hrBase(:,:,:,rep,:,:) = resp.(sess).hrBase(:,:,:,rep,perm,:);
%         end
%     end
%     respP{subjInd} = respP;
%     %Save
% %     save(fullfile(funPath,outDir,[subj '.mat']),'resP','-append')
% end
% 
% [respP, info] = reorgData(p,respP);
% 
% % disp('Permuting Labels: Done')