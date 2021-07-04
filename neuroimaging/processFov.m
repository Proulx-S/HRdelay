function filePath = processFov(p)

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
        outDir  = 'd'; if ~exist(fullfile(funPath,outDir),'dir'); mkdir(fullfile(funPath,outDir)); end
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp


%% Retinotopic feature selection
if p.featSel.fov.doIt && strcmp(p.featSel.fov.threshMethod,'empirical')
    filePath = fullfile(funPath,outDir,mfilename);
    forceFlag = 0;
    if ~exist([filePath '.mat'],'file') || forceFlag
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
        
        disp('-----------------------')
        disp('empiricalFov: Computing')
        [featSel_areaAndFov,cont,voxProp,pEmpirical,fAll,f] = empiricalFov(d,p);
        disp(['empiricalFov: saving to ' filePath])
        save(filePath,'featSel_areaAndFov','cont','voxProp','pEmpirical')
        for i = 1:length(fAll)
            saveas(fAll{i},[filePath '_' num2str(i) '.fig'])
            saveas(fAll{i},[filePath '_' num2str(i) '.jpg'])
        end
        saveas(f,[filePath '.fig'])
        saveas(f,[filePath '.jpg'])
        disp('empiricalFov: saved')
        disp('-----------------------')
    end
    disp(['empiricalFov: copying to ' fullfile(p.figOption.finalDir)])
    copyfile([filePath '.fig'],fullfile(p.figOption.finalDir))
end





