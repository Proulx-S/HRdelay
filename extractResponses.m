function extractResponses(p)
disp('------------------')
disp('extractResponses.m')
% For each session of each subject, extract BOLD response to the stimulus
% cycle by fitting a sinusoidal response model and a deconvolution model.

%% Define paths
dataIn = p.dataPath.V1;
dataOut = fullfile(p.dataPath.V1,'resp');
if ~exist(dataOut,'dir')
    mkdir(dataOut)
end

%% Run
for subjInd = 1:length(p.meta.subjList)
    % Load data
    subj = p.meta.subjList{subjInd};
    disp([subj ': loading time-series'])
    data = load(fullfile(dataIn,[subj '.mat']),'d','p');
    % Run GLMs
    for sessInd = 1:size(data.d.fun,2)
        disp([subj ': extracting responses from timeseries (sess' num2str(sessInd) '/' num2str(size(data.d.fun,2)) ')'])
        sess = ['sess' num2str(sessInd)];
        resp.(sess) = runGLMs(data.d.fun(1,sessInd),data.p,0);
        dDtrd.(sess) = rmfield(data.d.fun(sessInd),{'data' 'design' 'extraRegr'});
        dDtrd.(sess).data = resp.(sess).dataDtrd;
        resp.(sess) = rmfield(resp.(sess),'dataDtrd');
    end
    % Add info to res
    resp.sess1.voxProp = data.d.voxProp;
    resp.sess2.voxProp = data.d.voxProp;
    % Save
    disp([subj ': saving responses'])
    save(fullfile(dataOut,[subj '.mat']),'resp')
    disp([subj ': saved to ' fullfile(dataOut,[subj '.mat'])])
    clear resp
    disp([subj ': saving detrended timeseries'])
    save(fullfile(dataIn,[subj '_dDtrd.mat']),'dDtrd')
    disp([subj ': saved to ''' fullfile(dataIn,[subj '_dDtrd.mat']) ''''])
    clear dDtrd
end

