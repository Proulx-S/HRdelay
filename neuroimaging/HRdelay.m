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

plotDecoding_acc(res,figOption)
plotDecoding_auc(res,figOption)
plotDecoding_distT(res,figOption)



nPerm = 2^14;

svmSpace = 'cart';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartNoAmp';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartNoDelay';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartReal';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

accPerm = plotDecodingPerm_acc(res,figOption);
plotDecoding_acc(res,figOption,accPerm)

plotDecodingPerm_auc(res,figOption)

