p.wd = pwd;
p.anaID = '2021-08-11'; %an ID for your output folder. Potentially useful if you make tweaks

%% Dependencies
addpath(genpath(fullfile(p.wd,'utilities')));
addpath(genpath(fullfile(p.wd,'matlabFun')));

% Matlab proprietary toolboxes:
% -Image Processing
% -Statistics and Machine Learning
% -Mapping
% -Curve Fitting
% -Bioinformatics

% Matlab contributed toolboxes
matDependencyPath = fullfile(p.wd,'matlabFileEx'); if ~exist(matDependencyPath,'dir'); mkdir(matDependencyPath); end
%https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap
%https://www.mathworks.com/matlabcentral/fileexchange/41961-nanconv
addpath(genpath(fullfile(matDependencyPath)));


% gitDependencyPath = fullfile(p.paths.home,'Documents/GitHub/utilities');
%addpath(genpath(fullfile(gitDependencyPath,'bassFun/matlab')));
%addpath(genpath(fullfile(gitDependencyPath,'circstat-matlab')));
%addpath(genpath(fullfile(gitDependencyPath,'RAMBiNo')));
%addpath(genpath(fullfile(gitDependencyPath,'BrewerMap')));
% matDependencyPath = fullfile(p.paths.home,'OneDrive - McGill University\MATLAB');
% addpath(genpath(fullfile(matDependencyPath,'shadedErrorBar-master')));
%addpath(genpath(fullfile(matDependencyPath,'InterX')))
%addpath(genpath(fullfile(pwd,'matlabFun')));


p.dataPath.V1 = fullfile(p.wd,'data','V1');

%% Output paths
p.figOption.outDir = fullfile(p.wd,'figures',p.anaID); if ~exist(p.figOption.outDir,'dir'); mkdir(p.figOption.outDir); end
p.termOption.outDir = fullfile(p.wd,'termOutputs',p.anaID); if ~exist(p.termOption.outDir,'dir'); mkdir(p.termOption.outDir); end


%% Logging
p.figOption.verbose  = 0;
p.termOption.verbose = 0;
p.termOption.save = 1; % save command window outputs to text file
if p.termOption.save
    diaryON(p)
end


%% Input Data


% p.paths.repo.out = fullfile('McGill University/Farivar Lab - Dissertations/Sebastien/Manuscripts/aa - in preparation/SP_Neuroimage_HRdelay',p.anaID); if ~exist(p.paths.repo.out,'dir'); mkdir(p.paths.repo.out); end




%% Parameters for single-subject example
p.figOption.save     = 1; % save figures
p.figOption.subjInd  = 2;
p.figOption.sessInd  = 1;
p.figOption.condInd  = 1;
p.figOption.sliceInd = 7;



%% Permutation
p.perm.doIt = 1;
p.perm.n = 2^13;

%% Bootstrapping
p.boot.n = 2^13;



%% Meta data
p.meta.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
