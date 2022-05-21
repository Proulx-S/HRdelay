disp('--------------')
disp('downloadData.m')
%% Download intermediary data from Zenodo repository DOI:10.5281/zenodo.5183028 (https://doi.org/10.5281/zenodo.5183028)
forceDownload = 0;
if ~exist(p.dataPath.V1,'dir')
    mkdir(p.dataPath.V1)
end
outDir = fullfile(p.dataPath.V1,'ts');
if ~exist(outDir,'dir')
    mkdir(outDir)
end

for subjInd = 1:length(p.meta.subjList)
    curLink = ['https://zenodo.org/record/5183028/files/' p.meta.subjList{subjInd} '.mat?download=1'];
    curFile = fullfile(outDir,[p.meta.subjList{subjInd} '.mat']);
    if ~exist(curFile,'file') || forceDownload
        disp(['Downloading ' p.meta.subjList{subjInd} ' from https://doi.org/10.5281/zenodo.5183028'])
        websave(curFile,curLink,weboptions('timeout',15));
        disp(['Downloaded to ' curFile])
    else
        disp([p.meta.subjList{subjInd} ' already downloaded from https://doi.org/10.5281/zenodo.5183028'])
    end
end
