function processResponses(pMaster,figOption,verbose)
actuallyRun = 1;
if ~actuallyRun
    disp(['skipping ' mfilename])
    return
end
if ~exist('verbose','var')
    verbose = 1;
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
            inDir  = 'c';
            outDir = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Run
if actuallyRun
    for subjInd = 1:length(subjList)
        % Load data
        subj = subjList{subjInd};
        disp([subj ': loading'])
        load(fullfile(funPath,inDir,[subj '.mat']),'d','p')
        % Run GLMs
        for sessInd = 1:size(d.fun,2)
            disp([subj ': processing responses (sess' num2str(sessInd) '/' num2str(size(d.fun,2)) ')'])
            sess = ['sess' num2str(sessInd)];
            p.perm = pMaster.perm;
            res.(sess) = runGLMs(d.fun(1,sessInd),p,0);
            dDtrd.(sess) = rmfield(d.fun(sessInd),{'data' 'design' 'extraRegr'});
            dDtrd.(sess).data = res.(sess).dataDtrd;
            res.(sess) = rmfield(res.(sess),'dataDtrd');
        end
        % Add info to res
        res.sess1.voxProp = d.voxProp;
        res.sess2.voxProp = d.voxProp;
        db
        % Save
        disp([subj ': saving responses'])
        if ~exist(fullfile(funPath,outDir),'dir')
            mkdir(fullfile(funPath,outDir))
        end
        save(fullfile(funPath,outDir,[subj '.mat']),'res')
        clear res
        if ~pMaster.perm.doIt
            disp([subj ': saving detrended data'])
            if ~exist(fullfile(funPath,outDir),'dir')
                mkdir(fullfile(funPath,outDir))
            end
            save(fullfile(funPath,outDir,[subj '_dDtrd.mat']),'dDtrd')
            disp([subj ': saved to ''' fullfile(funPath,outDir) ''''])
        end
        clear dDtrd
    end
end

