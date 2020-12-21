clear all
addpath(genpath(fullfile(pwd,'matlabFun')));

runFit

fitType = 'mixed';
threshType = 'fdr';
veinPerc = 20;
figOption.save = 1;
figOption.subj = 1; % subjInd, +inf for all subj
% loadOption;
preprocAndShowMasks(fitType,threshType,veinPerc,figOption)

inspectSubjAndExclude(figOption)

runGroupAnalysis_sin(figOption)

runGroupAnalysis_hr(figOption)

clear res

svmSpace = 'cart_HT';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoAmp_HT';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoDelay_HT';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartReal_T';
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

