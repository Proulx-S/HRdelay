function processFov(p)
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end
if ~isfield(p,'perm')
    p.perm.doIt = 0;
end
if ~isfield(p,'termOption') || isempty(p.figOption)
    p.termOption.verbose = 1;
end


%% Define paths
subjList = p.meta.subjList;
repoPath = p.paths.repo.in;
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
dAll = cell(length(subjList),1);
for subjInd = 1:length(subjList)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if p.termOption.verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
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

%% Retinotopic feature selection
if p.featSel.fov.doIt && strcmp(p.featSel.fov.threshMethod,'empirical')
%     if p.perm.doIt
%         % Not actually using permuted data. But it does not acutally
%         % matters, because permuting does not change the result here. This
%         % step is base on the response delay irrespective of stimulus
%         % condition (delay computed after averaging across repetutions and
%         % stimulus conditions). Loading the precomputed empirical FOV
%         % therefore yields the same result as recomputing it, but save
%         % loads of time.
%         [featSel_fov,d,p] = empiricalFov(d,p,fullfile(funPath,outDir2));
%     else
        [featSel_fov,d,p] = empiricalFov(d,p,fullfile(funPath,outDir));
        % saves to empiricalFov.mat
%     end
end


%% Functional feature selection
featSel = cell(size(d));
f = cell(size(d));
disp('Feature Selection: processing')
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        p.subjInd = subjInd;
        p.sessInd = sessInd;
        disp(['subj' num2str(subjInd) '; sess' num2str(sessInd)])
        [featSel{subjInd,sessInd}] = getFeatSel(d{subjInd,sessInd},p,featSel_fov{subjInd,sessInd});
        featSel{subjInd,sessInd}.GLMs.random.sinDesign = d{subjInd,sessInd}.sinDesign;
        featSel{subjInd,sessInd}.GLMs.random.hrDesign = d{subjInd,sessInd}.hrDesign;
        featSel{subjInd,sessInd}.GLMs.random.infoDesign = d{subjInd,sessInd}.infoDesign;
        featSel{subjInd,sessInd}.GLMs.fixed.sinDesign = d{subjInd,sessInd}.featSel.F.act.design.full;
        
    end
end
% save to featSel.mat
disp('Feature Selection: saving')
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
save(fullfilename,'featSel')
disp(['Feature Selection: saved to ' fullfilename])

