disp('--------------')
disp('downloadData.m')
%% Download intermediary data (<600MB total) from Zenodo repository DOI:10.5281/zenodo.7058825 (https://doi.org/10.5281/zenodo.7058825)
forceDownload = 0;
if ~exist(p.dataPath.V1,'dir')
    mkdir(p.dataPath.V1)
end
outDir = fullfile(p.dataPath.V1,'ts');
if ~exist(outDir,'dir')
    mkdir(outDir)
end

for subjInd = 1:length(p.meta.subjList)
    curLink = ['https://zenodo.org/record/7058825/files/' p.meta.subjList{subjInd} '.mat?download=1'];
    curFile = fullfile(outDir,[p.meta.subjList{subjInd} '.mat']);
    if ~exist(curFile,'file') || forceDownload
        disp(['Downloading ' p.meta.subjList{subjInd} ' from https://doi.org/10.5281/zenodo.7058825'])
        websave(curFile,curLink,weboptions('timeout',15));
        disp(['Downloaded to ' curFile])
    else
        disp([p.meta.subjList{subjInd} ' already downloaded from https://doi.org/10.5281/zenodo.7058825'])
    end
end
