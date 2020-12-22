clear all
addpath(genpath(fullfile(pwd,'matlabFun')));

runFit

fitType = 'fixed';
threshType = 'p';
veinPerc = 20;
figOption.save = 1;
figOption.subj = inf; % subjInd, +inf for all subj
% loadOption;
preprocAndShowMasks(fitType,threshType,veinPerc,figOption)

inspectSubjAndExclude(figOption)

runGroupAnalysis_sin(figOption)

runGroupAnalysis_hr(figOption)

clear res

svmSpace = 'cart';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoAmp';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoDelay';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartReal';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

plotDecoding(res,figOption)
return

nPerm = 2^13;

svmSpace = 'cart_HT';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartNoAmp_HT';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartNoDelay_HT';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

plotDecodingPerm(res,figOption)

