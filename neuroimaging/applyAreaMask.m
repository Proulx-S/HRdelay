function applyAreaMask(figOption)
actuallyRun = 1;
if ~actuallyRun
    disp(['skipping ' mfilename])
    return
end
if exist('figOption','var')
    if ~isfield(figOption,'subj')
        figOption.subj = [];
    end
    if ~isfield(figOption,'save')
        figOption.save = 0;
    end
else
    figOption.subj = [];
    figOption.save = 0;
end

%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'b';
            outDir = 'c';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Define exclusion
exclude = 1;
exclusion.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
exclusion.subj = 2;
exclusion.sess = {1};
exclusion.run = {4};
exclusion.cond = {2};

%% Loop over subjects
for subjInd = 1:length(subjList)
    subj = subjList{subjInd};
    disp([subj ': loading'])
    load(fullfile(funPath,inDir,[subj '.mat']),'d','p')
    
    %% Store exclusion d
    if exclude
        subjIndX = ismember(exclusion.subj,subjInd);
        for sessInd = 1:size(d.fun,2)
            if any(subjIndX) && sessInd==exclusion.sess{subjIndX}
                d.fun(1,sessInd).excl = d.fun(1,sessInd).repLabel==exclusion.run{subjIndX};
            else
                d.fun(1,sessInd).excl = false(size(d.fun(1,sessInd).repLabel));
            end
        end
    end
    
    %% Run GLMs over the full cropped area for later figures
    if ismember(figOption.subj,subjInd)
        disp([subj ': fitting over the full crop area for later figure'])
        sessInd = 1;
        if any(ismember(figOption.subj,subjInd))
            res = runGLMs(d.fun(1,sessInd),p,0);
        end
        res = rmfield(res,'dataDtrd');
        % save
        if ~exist(fullfile(funPath,outDir),'dir')
            mkdir(fullfile(funPath,outDir))
        end
        save(fullfile(funPath,outDir,[subj '_fullFit.mat']),'res')
        clear res
    end
    
    %% Apply area mask and vectorize
    v1mask = false(size(d.fun(1,1).data{1},[1 2 3]));
    v1mask(:) = p.masks.roiMasks.v1(p.masks.cropMask);
    for sessInd = 1:size(d.fun,2)
        for runInd = 1:size(d.fun(1,sessInd).data,1)
            tmp = permute(d.fun(1,sessInd).data{runInd},[4 1 2 3]);
            d.fun(1,sessInd).data{runInd} = permute(tmp(:,v1mask),[2 3 4 1]); clear tmp
        end
    end
    %% Save
    if ~exist(fullfile(funPath,outDir),'dir')
        mkdir(fullfile(funPath,outDir))
    end
    disp([subj ': saving masked data'])
    save(fullfile(funPath,outDir,[subj '.mat']),'d','p')
    disp([subj ': saved to ''' fullfile(funPath,outDir) ''''])
    clear d p
end
