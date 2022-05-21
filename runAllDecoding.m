function [resBS,resBShr,resWS,f,info,fullfilenameDecoding] = runAllDecoding(p,permFlag,respP,featSelP)
if ~exist('permFlag','var')
    permFlag = 0;
end
if ~exist('respP','var') || isempty(respP)
    respP = [];
end
if ~exist('featSelP','var') || isempty(featSelP)
    featSelP = [];
end

if permFlag
    p.figOption.verbose = 0;
    p.termOption.verbose = 0;
end



%% Define paths
paths.subjList = p.meta.subjList;
paths.funPath = p.dataPath.V1;
paths.inDir = p.dataPath.V1;
paths.inDir2 = p.dataPath.V1;
paths.outDir = p.dataPath.V1;


condPairList = p.svm.condPairList;
respFeatList = p.svm.respFeatList;


figOption_verbose = 1;
resBS = cell(length(condPairList),length(respFeatList));
resBShr = cell(length(condPairList),length(respFeatList));
resWS = cell(length(condPairList),length(respFeatList));
f = cell(length(condPairList),length(respFeatList));
p.complexSpace = p.svm.complexSpace;
for respFeatInd = 1:length(respFeatList)
    p.chanSpace = respFeatList{respFeatInd};
    for condPairInd = 1:length(condPairList)
        if condPairInd~=1
            p.figOption.verbose = 0;
        else
            p.figOption.verbose = figOption_verbose;
        end
        p.condPair = condPairList{condPairInd};
        [resBS{condPairInd,respFeatInd},resBShr{condPairInd,respFeatInd},resWS{condPairInd,respFeatInd},f{condPairInd,respFeatInd}] = runDecoding(p,permFlag,paths,respP,featSelP);
    end
end
p.figOption.verbose = figOption_verbose;
info.condPairList = condPairList';
info.respFeatList = respFeatList;
info.info = 'condPair x respFeat';

if ~permFlag
    %% Save data
    fullpath = fullfile(p.wd,'results',p.anaID);
    if ~exist(fullpath,'dir'); mkdir(fullpath); end
    fullfilenameDecoding = fullfile(fullpath,'decoding');
    save(fullfilenameDecoding,'resBS','p','info');
    fullfilename = fullfile(fullpath,'channels');
    save(fullfilename,'resBShr','p','info');
    disp(['decoding results saved to ' fullfilename])
end

