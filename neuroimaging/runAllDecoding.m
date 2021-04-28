function [resBS,resBShr,resWS,f,info] = runAllDecoding(p,condPairList,respFeatList,figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
if ~exist('condPairList','var') || isempty(condPairList)
    condPairList = {'grat1VSgrat2' 'grat1VSplaid' 'grat2VSplaid'};
end
if ~exist('respFeatList','var') || isempty(respFeatList)
    respFeatList = {'cart' 'cartNoDelay' 'cartNoAmp'};
end



figOption_verbose = p.figOption.verbose;


resBS = cell(length(condPairList),length(respFeatList));
resBShr = cell(length(condPairList),length(respFeatList));
resWS = cell(length(condPairList),length(respFeatList));
f = cell(length(condPairList),length(respFeatList));
for respFeatInd = 1:length(respFeatList)
    p.svmSpace = respFeatList{respFeatInd};
    for condPairInd = 1:length(condPairList)
        if condPairInd~=1
            p.figOption.verbose = 0;
        else
            p.figOption.verbose = figOption_verbose;
        end
        p.condPair = condPairList{condPairInd};
        [resBS{condPairInd,respFeatInd},resBShr{condPairInd,respFeatInd},resWS{condPairInd,respFeatInd},f{condPairInd,respFeatInd}] = runDecoding(p,verbose);
    end
end
info.condPairList = condPairList';
info.respFeatList = respFeatList;
info.info = 'condPair x respFeat';


