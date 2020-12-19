clear all
addpath(genpath(fullfile(pwd,'matlabFun')));

runFit

fitType = 'fixed';
threshType = 'p';
veinPerc = 20;
figOption.save = 1;
figOption.subj = inf; % subjInd, +inf for all subj
% loadOption;
defineAndShowMasks(fitType,threshType,veinPerc,figOption)

applyFeatSelAndClean(figOption)

runGroupAnalysis_sin(figOption)

runGroupAnalysis_hr(figOption)

clear res

svmSpace = 'cart_HT';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'cartNoAmp_HT';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

svmSpace = 'polMag_T';
resTmp = runDecoding(svmSpace);
res.(svmSpace) = resTmp;

plotDecoding(res,figOption)
return

nPerm = 1000;

svmSpace = 'cart_HT';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'cartNoAmp_HT';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

svmSpace = 'polMag_T';
res.(svmSpace) = runDecoding(res.(svmSpace),nPerm);

plotDecodingPerm(res,figOption)

% %%javascript
% IPython.notebook.save_notebook()
% 
% %%shell
% jupyter nbconvert --to html HRdelay.ipynb
