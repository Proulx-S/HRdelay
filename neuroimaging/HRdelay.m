clear all
% Dependencies
addpath(genpath(fullfile(pwd,'matlabFun')));
verbose = 0;

% Between-session feature selection parameters
% activated voxels
featSel_bSess.activation.doIt = 1;
featSel_bSess.activation.fitType = 'fixed';
featSel_bSess.activation.threshType = 'p';
featSel_bSess.activation.threshVal = 0.05;
% vein voxels
featSel_bSess.vein.doIt = 1;
featSel_bSess.vein.source = 'fullModelResid';% 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
featSel_bSess.vein.percentile = 20;
% discriminant voxels
featSel_bSess.discrim.doIt = 1;
% Display parameters
figOption.save = 1;
figOption.subj = inf; % subjInd, +inf for all subj



runFit(verbose)

preprocAndShowMasks(featSel_bSess,figOption,verbose)

inspectSubjAndExclude(figOption,verbose)

runGroupAnalysis_sin(figOption,verbose)

runGroupAnalysis_hr(figOption,verbose)

runAllDecoding(figOption,verbose)

