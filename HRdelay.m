%% Info
% This script reproduces figures from S. Proulx et al., Cerebral Cortex, in
% preperation.
figure('MenuBar','none','ToolBar','none');
clear all
close all

%% Initiation
initAnalysis; % where to define paths, analysis parameters and some other general configurations

%% Get data (preprocessed BOLD volumes, only V1 ROI voxels, one .mat file per subject, total<600MB)
downloadData;

%% Extract sinusoidal and model-free responses from timeseries
extractResponses(p);

%% Define the retinotopic representation of the stimulus field of view (data-driven with priors from a probabilistic retinotopic atlas)
doWhat = 'download';
% 'run'             -> run it yourself, but will not run when locally saved data is available (several minutes; uses some more Matlab proprietary toolboxes, see initAnalysis.m)
% 'run_forced'      -> same as above, but forces to rerun instead of loading locally saved data
% 'download'        -> download precomputed data from repository, but will not download when locally saved data is available
% 'download_forced' -> same as above, but forces to redownload instead of loading locally saved data
% 'run_forced_butSkipFlattenEccDist' -> same as run_forced, but skips forced run of flattenEccDist (very long step)
processFov(p,doWhat);

% % Figure 3A histogram inset
% open('/Users/sebastienproulx/HRdelay/figures/myAnalysis/Fig3Ahist.fig')
% f = gcf;
% ax = gca;
% ax.Children(3).FaceColor = 'k';
% ax.Children(3).EdgeColor = 'none';
% f.Color = 'none';
% ax.Color = 'none';
% ax.Box = 'off';
% ax.YAxis.Visible = 'off';
% ax.XTick = [];
% % saveas(f, 'Fig3Ahist.svg')

% % Figure SuppFig1B histograms
% open('/Users/sebastienproulx/HRdelay/figures/myAnalysis/SuppFig1hist.fig')
% f = gcf;
% ax = f.Children;
% hHist = findobj([ax.Children],'Type','Histogram');
% set(hHist,'FaceColor','k','EdgeColor','none');
% f.Color = 'none';
% set(ax,'Color','none');
% set(ax,'XTick',-pi/2:pi/2:(pi+pi/2));
% set(ax,'XTickLabel',-3:3:9);
% % saveas(f, 'SuppFig1hist.svg')



return

%% Feature selection
processFeatSel(p);

%% Decoding
runAllDecoding(p);

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
