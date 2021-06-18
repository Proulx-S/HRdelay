function [featSelP,featSel_fov] = processFeatSel(p,resPall,featSel_fov)
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end
if ~isfield(p,'termOption') || isempty(p.figOption)
    p.termOption.verbose = 1;
end
if ~exist('resPall','var')
    permFlag = 0;
else
    permFlag = 1;
end
if permFlag
    p.figOption.verbose = 0;
    p.termOption.verbose = 0;
    verbose = 0;
else
    verbose = 1;
end

%% Define paths
subjList = p.meta.subjList;
repoPath = p.paths.repo.in;
    funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
        inDir  = 'd';
        outDir  = ['e_' p.anaID];
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:length(subjList)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if p.termOption.verbose; disp(['loading: ' curFile]); end
    if ~permFlag
        load(curFile,'res');
    else
        res = resPall{subjInd}; resPall{subjInd} ={};
%         load(curFile,'resP');
%         res = resP; clear resP
    end
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
if ~exist('featSel_fov','var') || isempty(featSel_fov)
    if p.featSel.fov.doIt && strcmp(p.featSel.fov.threshMethod,'empirical')
        filePath = fullfile(funPath,inDir,'processFov');
        load(filePath,'cont','voxProp','pEmpirical','featSel_areaAndFov')
        
        %Pack cont into featSel_areaAndFov and voxProp into d
        for subjInd = 1:size(d,1)
            for sessInd = 1:size(d,2)
                featSel_areaAndFov{subjInd,sessInd}.cont.L = cont.L{subjInd,sessInd};
                featSel_areaAndFov{subjInd,sessInd}.cont.R = cont.R{subjInd,sessInd};
                d{subjInd,sessInd}.voxProp.L = voxProp{subjInd,sessInd}.L;
                d{subjInd,sessInd}.voxProp.R = voxProp{subjInd,sessInd}.R;
            end
        end
        %Pack pEmirical into p
        p.featSel.fov.empricalFov = pEmpirical;
        
        featSel_fov = featSel_areaAndFov; clear featSel_areaAndFov
    end
end

%% Functional feature selection
featSel = cell(size(d));
f = cell(size(d));
if verbose
    disp('Feature Selection: processing')
end
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        if verbose
            disp(['subj' num2str(subjInd) '; sess' num2str(sessInd)])
        end
        [featSel{subjInd,sessInd}] = getFeatSel(d{subjInd,sessInd},p,featSel_fov{subjInd,sessInd});
        featSel{subjInd,sessInd}.GLMs.random.sinDesign = d{subjInd,sessInd}.sinDesign;
        featSel{subjInd,sessInd}.GLMs.random.hrDesign = d{subjInd,sessInd}.hrDesign;
        featSel{subjInd,sessInd}.GLMs.random.infoDesign = d{subjInd,sessInd}.infoDesign;
        featSel{subjInd,sessInd}.GLMs.fixed.sinDesign = d{subjInd,sessInd}.featSel.F.act.design.full;
    end
end
% save to featSel.mat
if verbose
    disp('Feature Selection: saving')
end
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
if ~permFlag
    save(fullfilename,'featSel')
    featSelP = [];
else
    featSelP = featSel; % clear featSel;
%     save(fullfilename,'featSelP','-append')
end
if verbose
    disp(['Feature Selection: saved to ' fullfilename])
end

