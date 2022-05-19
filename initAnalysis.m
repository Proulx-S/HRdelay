disp('--------------')
disp('initAnalysis.m')

p.wd = pwd;
p.anaID = 'anaTmp'; %an ID for your output folder. Potentially useful if you make tweaks

%% Dependencies
addpath(genpath(fullfile(p.wd,'fun')));

% Matlab proprietary toolboxes (may not all be needed):
% -Statistics and Machine Learning

% -Image Processing
% -Mapping
% -Curve Fitting
% -Bioinformatics

% Matlab community toolboxes
matDependencyPath = fullfile(p.wd,'matlabFileExchange'); if ~exist(matDependencyPath,'dir'); mkdir(matDependencyPath); end
matDependencyNameList = {...
    'https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap'...
    };
%https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap
%https://www.mathworks.com/matlabcentral/fileexchange/41961-nanconv
%https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
matDependencyUrlList = {...
    'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/25536/versions/1/download/zip'...
    };
for i = 1:length(matDependencyNameList)
    curLink = matDependencyUrlList{i};
    curZip = fullfile(matDependencyPath,'x.zip');
    websave(curZip,curLink,weboptions('timeout',15));
    [~,b,~] = fileparts(matDependencyNameList{i});
    curFolder = fullfile(matDependencyPath,b);
    unzip(curZip,curFolder)
    delete(curZip)
    addpath(genpath(curFolder))
end

% Github repository (make sure git is installed and available in your system path: https://git-scm.com/downloads)
gitDependencyPath = fullfile(p.wd,'..');
% libsvm
curRepo = 'libsvm';
curRepoURL = 'https://github.com/cjlin1/libsvm.git';
curRepoPath = fullfile(gitDependencyPath,curRepo);
disp(curRepoURL)
if exist(fullfile(curRepoPath,'.git'),'dir')
    disp(['already in ' curRepoPath])
else
    eval(['!git clone ' curRepoURL ' ' curRepoPath])
    if ispc
        addpath(genpath(fullfile(curRepoPath,'windows')));
    else
        error('you need to fix paths to libsvm githu repo for non-windows machines')
    end
end
% shadedErrorBar
curRepo = 'shadedErrorBar';
curRepoURL = 'https://github.com/raacampbell/shadedErrorBar.git';
curRepoPath = fullfile(gitDependencyPath,curRepo);
disp(curRepoURL)
if exist(fullfile(curRepoPath,'.git'),'dir')
    disp(['already in ' curRepoPath])
else
    eval(['!git clone ' curRepoURL ' ' curRepoPath])
    addpath(genpath(curRepoPath));
end

%% Parameters for feature selection
% Within visual field region of the stimulus
p.featSel.fov.doIt = 1;
p.featSel.fov.threshMethod = 'empirical'; % 'empirical' 'ecc'
p.featSel.fov.areaLabel = 'v1';
p.featSel.fov.threshVal = [0.75 7]; % threshMethod='ecc'
p.featSel.fov.percentile = 20; % threshMethod='ecc'
% Activated voxels (random-effect)
p.featSel.respVecSig.doIt = 1;
p.featSel.respVecSig.threshMethod = 'p'; % '%ile' 'p' 'fdr'
p.featSel.respVecSig.threshVal = 0.05; % threshMethod='p' or 'fdr'
p.featSel.respVecSig.percentile = 20; % threshMethod='%ile'
% Less likely-to-be-vein voxels
p.featSel.vein.doIt = 1;
p.featSel.vein.threshMethod = '%ile'; % '%ile'
p.featSel.vein.threshVal = 0.01; % not used
p.featSel.vein.percentile = 20; % threshMethod='%ile'
% Most discriminant voxels
p.featSel.respVecDiff.doIt = 1;
p.featSel.respVecDiff.threshMethod = '%ile'; % '%ile' 'p' 'fdr'
p.featSel.respVecDiff.threshVal = 0.5; % threshMethod='p' or 'fdr'
p.featSel.respVecDiff.percentile = 20; % threshMethod='%ile'
% Feature Combination
p.featSel.global.doIt = 1;
p.featSel.global.method = 'allData';
% 'allData'-> featSel uses all three conditions, irrespective of the condition pairs to be decoded
% 'custom1'-> featSel of active voxels uses all three conditions but featSel of discriminant voxels uses only the conditions to be decoded
% 'custom2'-> featSel of active and most discriminant voxels uses only the conditions to be decoded


%% SVM parameters
p.svm.kernel.type = 'lin';
p.svm.complexSpace = 'bouboulisDeg1'; % 'bouboulisDeg1' 'bouboulisDeg2'
p.svm.doWithin = 0;

p.svm.condPairList = {'grat1VSgrat2' 'grat1VSplaid' 'grat2VSplaid'};
p.svm.respFeatList = {'delay' 'cartNoDelay' 'cart'};


%% Meta data
p.meta.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';


%% Paths
% p.figPath = fullfile(p.wd,'figures'); if ~exist(p.figPath,'dir'); mkdir(p.figPath); end
p.figOption.outDir = fullfile(p.wd,'figures',p.anaID); if ~exist(p.figOption.outDir,'dir'); mkdir(p.figOption.outDir); end
p.termOption.outDir = fullfile(p.wd,'terminalOutputs',p.anaID); if ~exist(p.termOption.outDir,'dir'); mkdir(p.termOption.outDir); end
p.dataPath.V1 = fullfile(p.wd,'data','V1');


%% Logging
% p.figOption.verbose  = 2; % 0, 1 or 2
% p.termOption.verbose = 2;

% p.termOption.save = 1; % save command window outputs to text file
% if p.termOption.save
    diaryON(p)
% end


%% Parameters for single-subject example
% p.figOption.save     = 1; % save figures
p.figOption.subjInd  = 2;
p.figOption.sessInd  = 1;
p.figOption.condInd  = 1;
p.figOption.sliceInd = 7;


%% Permutations
p.perm.doIt = 1;
p.perm.n = 2^13;


%% Bootstrapping
p.boot.n = 2^13;
