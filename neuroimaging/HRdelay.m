clear all
close all
%% Dependencies
gitDependencyPath = '/Users/sebastienproulx/Documents/GitHub/utilities';
matDependencyPath = '/Users/sebastienproulx/Dropbox/MATLAB';
addpath(genpath(fullfile(gitDependencyPath,'circstat-matlab')));
addpath(genpath(fullfile(gitDependencyPath,'RAMBiNo')));
addpath(genpath(fullfile(pwd,'matlabFun')));
% rmpath(genpath(fullfile(matDependencyPath,'HotellingT2')));
verbose = 1; % prints more info

%% Feature selection parameters
% Within visual field region of the stimulus
p.featSel.fov.doIt = 1;
p.featSel.fov.threshMethod = 'empirical'; % 'empirical' 'ecc'
p.featSel.fov.areaLabel = 'v1';
p.featSel.fov.threshVal = [0.75 7]; % threshMethod='ecc'
p.featSel.fov.percentile = 20; % threshMethod='ecc'
% Most activated voxels
p.featSel.act.doIt = 1;
p.featSel.act.threshMethod = 'fdr'; % '%ile' 'p' 'fdr'
p.featSel.act.threshVal = 0.05; % threshMethod='p' or 'fdr'
p.featSel.act.percentile = 20; % threshMethod='%ile'
% Most significant response vectors
p.featSel.respVecSig.doIt = 1;
p.featSel.respVecSig.threshMethod = 'fdr'; % '%ile' 'p' 'fdr'
p.featSel.respVecSig.threshVal = 0.05; % threshMethod='p' or 'fdr'
p.featSel.respVecSig.percentile = 20; % threshMethod='%ile'
% Less likely-to-be-vein voxels
p.featSel.vein.doIt = 1;
p.featSel.vein.threshMethod = '%ile'; % '%ile'
p.featSel.vein.threshVal = 0.01; % not used
p.featSel.vein.percentile = 20; % threshMethod='%ile'
% Most discriminant voxels
p.featSel.respVecDiff.doIt = 1;
p.featSel.respVecDiff.threshMethod = '%ile'; % '%ile' 'p' 'fdr'
p.featSel.respVecDiff.threshVal = 0.5; % threshMethod='p' or 'fdr'
p.featSel.respVecDiff.percentile = 20; % threshMethod='%ile'
% Feature Combination
p.featSel.global.doIt = 1;
p.featSel.global.method = 'custom2';
% 'allData'-> featSel uses all three conditions, irrespective of the condition pairs to be decoded
% 'custom1'-> featSel of active voxels uses all three conditions but featSel of discriminant voxels uses only the conditions to be decoded
% 'custom2'-> featSel of active and most discriminant voxels uses only the conditions to be decoded

%% Normalization parameters
% p.norm.doCartSpaceScale = 1;
%% Display parameters
figOption.save = 0; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)
p.figOption.verbose  = 1;
p.figOption.subjInd  = figOption.subj;
p.figOption.sessInd  = 1;
p.figOption.sliceInd = 10;
if 0
    importData(verbose)
    applyAreaMask(figOption)
    processResponses(figOption,verbose)
%     processWaveletResponses(figOption,verbose)
end
processFeatSel(p,verbose)
visualizeFeatSel(p)
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
