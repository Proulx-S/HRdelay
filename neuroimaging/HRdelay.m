clear all
%% Dependencies
dependencyPath = '/Users/sebastienproulx/Documents/GitHub/utilities';
addpath(genpath(fullfile(dependencyPath,'circstat-matlab')));
addpath(genpath(fullfile(pwd,'matlabFun')));
verbose = 1; % prints more info

%% Feature selection parameters
% Less likely-to-be-vein voxels
p.featSel.vein.doIt = 1;
p.featSel.vein.threshMethod = '%ile'; % '%ile' 'p' 'fdr'
p.featSel.vein.threshVal = 0.01;
p.featSel.vein.percentile = 20;
% Most activated voxels
p.featSel.act.doIt = 1;
p.featSel.act.threshMethod = 'fdr'; % '%ile' 'p' 'fdr'
p.featSel.act.threshVal = 0.05;
p.featSel.act.percentile = 20;
% Most significant response vectors
p.featSel.respVecSig.doIt = 1;
p.featSel.respVecSig.threshMethod = 'fdr'; % '%ile' 'p' 'fdr'
p.featSel.respVecSig.threshVal = 0.05;
p.featSel.respVecSig.percentile = 20;
% Most discriminant voxels
p.featSel.discrim.doIt = 1;
p.featSel.discrim.threshMethod = '%ile'; % '%ile' 'p' 'fdr'
p.featSel.discrim.threshVal = 0.03;
p.featSel.discrim.percentile = 20;
% Feature Combination
p.featSel.global.doIt = 1;
p.featSel.global.method = 'all';
% p.featSel.global.percentile = 30;

%% Normalization parameters
p.norm.doCartSpaceScale = 1;
%% Display parameters
figOption.save = 0; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)
p.figOption.verbose = 1;
p.figOption.subjInd = figOption.subj;
p.figOption.sessInd = 1;

if 0
    importData(verbose)
    applyAreaMask(figOption)
    processResponses(figOption,verbose)
    processWaveletResponses(figOption,verbose)
    processFeatSel(p,verbose)
end
visualizeFeatSel(p)
return
[resBS,resWS] = runAllDecoding(p,figOption,verbose);
groupAna(p,figOption,verbose)



return

% runFit(verbose,figOption)
runFit2(verbose,figOption)
runWave2(verbose,figOption)

preprocAndShowMasks(p,figOption,verbose)
inspectSubjAndExclude(figOption,verbose)



runGroupAnalysis_sin(figOption,verbose)
runGroupAnalysis_hr(figOption,verbose)

verbose = 1; % prints more info

return
runAllDecodingPerm(res,figOption,verbose);
