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
p.featSel.vein.threshVal = 0.001;
p.featSel.vein.percentile = 80;
% Most activated voxels
p.featSel.act.doIt = 1;
p.featSel.act.threshMethod = 'fdr';
p.featSel.act.threshVal = 0.001;
p.featSel.act.percentile = 20;
% Most significant response vectors
p.featSel.respVectSig.doIt = 1;
p.featSel.respVectSig.threshMethod = 'fdr';
p.featSel.respVectSig.threshVal = 0.001;
p.featSel.respVectSig.percentile = 20;
% Most discriminant voxels
p.featSel.discrim.doIt = 1;
p.featSel.discrim.threshMethod = '%ile';
p.featSel.discrim.threshVal = 0.001;
p.featSel.discrim.percentile = 20;
%% Normalization parameters
p.norm.doCartSpaceScale = 1;
%% Display parameters
figOption.save = 0; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)

if 0
    importData(verbose)
    applyAreaMask(figOption)
    processResponses(figOption,verbose)
    processWaveletResponses(figOption,verbose)
end

[resBS,resWS] = runAllDecoding(p,figOption,verbose);
return
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
