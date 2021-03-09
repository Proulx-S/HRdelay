function processFeatSel(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,2)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
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


%% Feature selection
featSel = cell(size(d));
featSel2 = cell(size(d));
tic
parfor i = 1:numel(d)
        [subjInd,sessInd] = ind2sub(size(d),i);
        [featSel{i},featSel2{i}] = getFeatSel(d{subjInd,sessInd},p,subjInd,sessInd);
end
toc
if p.figOption.verbose>=1
    subjInd = p.figOption.subjInd;
    sessInd = p.figOption.sessInd;
    x = featSel2{subjInd,sessInd}.percAll(:,1);
    y = featSel2{subjInd,sessInd}.percAll(:,2);
    z = featSel2{subjInd,sessInd}.percAll(:,3);
    indIn = featSel2{subjInd,sessInd}.indIn;
    featLabel = featSel2{subjInd,sessInd}.featLabel;
    
    figure('WindowStyle','docked');
    scatter3(x(indIn),y(indIn),z(indIn),'k.'); hold on
    scatter3(x(~indIn),y(~indIn),z(~indIn),'r.');
    % ax = gca;
    % ax.CameraPosition = camPos;
    xlabel(['' featLabel{1} ''])
    ylabel(['' featLabel{2} ''])
    zlabel(['' featLabel{3} ''])
    ax = gca;
    ax.XAxis.Scale = 'log';
    ax.YAxis.Scale = 'log';
    ax.ZAxis.Scale = 'log';
%     xlabel(['zscore(-log(' info1{1+1} '))'])
%     ylabel(['zscore(log(' info1{1+2} '))'])
%     zlabel(['zscore(log(' info1{1+3} '))'])
    legend({'include' 'exclude'})
end


%% Save
disp('saving feature selection')
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
save(fullfilename,'featSel','featSel2')
disp(['saved to: ' fullfilename])
