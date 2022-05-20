%% Info
% This script reproduces all figures from S. Proulx et al., Neuroimage, in
% preperation.
clear all
close all

%% Initiation
initAnalysis; % where to define paths, analysis parameters and some other general configurations

%% Get data (preprocessed BOLD volumes, only V1 ROI voxels, one .mat file per subject)
downloadData;

%% Extract responses from timeseries
extractResponses(p);

%% Define the retinotopic representation of the stimulus field of view (data-driven with priors from a probabilistic retinotopic atlas)
doWhat = 'download';
% 'run' -> run it yourself (several minutes; and will likely use some of Matlab's proprietary toolboxes)
% 'download' -> download precomputed data from repository
% 'run_forced' -> same as above, but forces to rerun instead of loading
% locally saved data
% 'download_forced' -> same as above, but forces to redownload instead of loading
% locally saved data
processFov(p,doWhat);

%% Feature selection
processFeatSel(p);

%% Decoding
runAllDecoding(p);
% [resBS,resBShr,resWS,f,info,decodingOut] = runAllDecoding(p);

%% Visualize and print stats on BOLD responses
visualizeResponses(p)

%% Visualize and print stats on decoding
plotAllDecoding(p);
statsAllDecoding(p);

%% Run permutation test (will take long time)
runPermDecoding(p)
plotAllDecoding(p);
return




%% Old code that will not run
chan = processChanHr(p,resBShr,info);
f = plotChanHr(p,chan);
statsChanHr(p,chan);
groupAna(p,figOption,verbose)
% runFit(verbose,figOption)
runFit2(verbose,figOption)
runWave2(verbose,figOption)
preprocAndShowMasks(p,figOption,verbose)
inspectSubjAndExclude(figOption,verbose)
runGroupAnalysis_sin(figOption,verbose)
runGroupAnalysis_hr(figOption,verbose)
