clear all
%% Dependencies
addpath(genpath(fullfile(pwd,'matlabFun')));
verbose = 1; % prints more info

%% Between-session feature selection parameters
% Most activated voxels
p.act.doIt = 1;
p.act.threshVal = 0.00000001;
p.act.percentile = 50;
% Less likely-to-be-vein voxels
p.vein.doIt = 1;
p.vein.percentile = 0;
% Most discriminant voxels
p.discrim.doIt = 1;
p.discrim.percentile = 0;
%% Display parameters
figOption.save = 0; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)

if 0
    importData(verbose)
    applyAreaMask(figOption)
    processResponses(figOption,verbose)
end
processWaveletResponses(figOption,verbose)
return

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
