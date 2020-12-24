clear all
addpath(genpath(fullfile(pwd,'matlabFun')));

runFit(0)


featSel_bSess.activation.doIt = 1;
featSel_bSess.activation.fitType = 'fixed';
featSel_bSess.activation.threshType = 'p';
featSel_bSess.activation.threshVal = 0.05;
featSel_bSess.vein.doIt = 1;
featSel_bSess.vein.source = 'fullModelResid';% 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
featSel_bSess.vein.percentile = 20;
featSel_bSess.discrim.doIt = 1;
% featSel_bSess.discrim.spaceList = {'cart'};

figOption.save = 0;
figOption.subj = 1; % subjInd, +inf for all subj
% loadOption;


preprocAndShowMasks(featSel_bSess,figOption,0)

inspectSubjAndExclude(figOption,0)

runGroupAnalysis_sin(figOption,0)

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

[P,H,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right'); P
disp('cart VS cartReal:')
disp(['signed rank = ' num2str(STATS.signedrank)])
disp(['one-sided p = ' num2str(P)])

plotDecoding_acc(res,figOption)


nPerm = 2^14; % will not run if already run, and instead just load previous results

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

