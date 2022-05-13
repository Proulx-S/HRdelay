function extractResponses(p)
actuallyRun = 1;
% if ~exist('verbose','var')
%     verbose = 1;
% end
% if exist('figOption','var')
%     if ~isfield(figOption,'subj')
%         figOption.subj = [];
%     end
%     if ~isfield(figOption,'save')
%         figOption.save = 0;
%     end
% else
%     figOption.subj = [];
%     figOption.save = 0;
% end


% %% Define paths
% subjList = p.meta.subjList;
% repoPath = p.paths.repo.in;
%         funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
%             inDir  = 'c';
% %             outDir = ['d_' p.anaID];
%             outDir = 'd';
% %make sure everything is forward slash for mac, linux pc compatibility
% for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
%     eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
% end
% clear tmp

%% Run
subjList = p.info.subjList;
funPath = p.paths.data;
inDir = '';
for subjInd = 1:length(subjList)
    % Load data
    subj = subjList{subjInd};
    disp([subj ': loading'])
    data = load(fullfile(funPath,inDir,[subj '.mat']),'d','p');
    % Run GLMs
    for sessInd = 1:size(data.d.fun,2)
        disp([subj ': processing responses (sess' num2str(sessInd) '/' num2str(size(data.d.fun,2)) ')'])
        sess = ['sess' num2str(sessInd)];
        %             p.perm = pMaster.perm;
        res.(sess) = runGLMs(data.d.fun(1,sessInd),data.p,0);
        dDtrd.(sess) = rmfield(data.d.fun(sessInd),{'data' 'design' 'extraRegr'});
        dDtrd.(sess).data = res.(sess).dataDtrd;
        res.(sess) = rmfield(res.(sess),'dataDtrd');
    end
    % Add info to res
    res.sess1.voxProp = data.d.voxProp;
    res.sess2.voxProp = data.d.voxProp;
    
    % Save
    outDir = '';
    filename = fullfile(funPath,outDir,[subj '_res.mat']);
    disp([subj ': saving responses to ' filename])
    save(filename,'res')
    clear res

    filename = fullfile(funPath,outDir,[subj '_dDtrd.mat']);
    disp([subj ': saving detrended data to ' filename])
    save(filename,'dDtrd')
    clear dDtrd
end

