clear all
close all

%% Display parameters
figOption.save = 1; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)
p.figOption.subjInd  = figOption.subj;
p.figOption.sessInd  = 1;
p.figOption.sliceInd = 7;
p.figOption.verbose  = 2;
p.figOption.save  = 2;
p.termOption.verbose = 1;
p.termOption.save = 0;

%% Open diary
if p.termOption.verbose && p.termOption.save
    eval(['diary ' fullfile(pwd,[mfilename '-' datestr(now,'YYYYmmDD-hh_MM_ss') '.log'])])
    disp(fullfile(pwd,[mfilename '.mat']))
    disp(datestr(now))
end

%% Dependencies
gitDependencyPath = '/Users/sebastienproulx/Documents/GitHub/utilities';
matDependencyPath = '/Users/sebastienproulx/Dropbox/MATLAB';
addpath(genpath(fullfile(gitDependencyPath,'circstat-matlab')));
addpath(genpath(fullfile(gitDependencyPath,'RAMBiNo')));
addpath(genpath(fullfile(pwd,'matlabFun')));
% rmpath(genpath(fullfile(matDependencyPath,'HotellingT2')));
verbose = 1; % prints more info

%% Meta data
p.meta.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';


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

%% SVM channel parameters
p.svm.condPairList = {'grat1VSgrat2' 'grat1VSplaid' 'grat2VSplaid'};
p.svm.respFeatList = {'cart' 'cartNoDelay' 'cartNoAmp'};

%% Display parameters
figOption.save = 1; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)
p.figOption.subjInd  = figOption.subj;
p.figOption.sessInd  = 1;
p.figOption.sliceInd = 10;
p.figOption.verbose  = 2;
p.figOption.save  = 2;
p.termOption.verbose = 1;
p.termOption.save = 1;

if 0
    importData(verbose)
    applyAreaMask(figOption)
    processResponses(p,figOption,verbose)
    %     processWaveletResponses(figOption,verbose)
end
if 1
    processFeatSel(p)
end
if 0
    visualizeFeatSel(p)
end
if 1
    [resBS,resBShr,resWS,f,info] = runAllDecoding(p,verbose);
    plotAllDecoding(p,resBS,info)
    statsAllDecoding(p,resBS,info)
end
return
if 0
    chan = processChanHr(p,resBShr,info);
    f = plotChanHr(p,chan);
    statsChanHr(p,chan);
end
if 1
    pPerm = p;
    pPerm.perm.doIt = 1;
    pPerm.perm.n = 5;
    pPerm.figOption.verbose = 0;
    pPerm.svm.condPairList = p.svm.condPairList(1);
    resPerm.perfMetric = nan(length(pPerm.meta.subjList),length(pPerm.svm.condPairList),length(pPerm.svm.respFeatList),pPerm.perm.n);
    for permInd = 1:pPerm.perm.n
        tic
        processResponses(pPerm,figOption,verbose)
        toc
        processFeatSel(pPerm)
        toc
        [resBS,~,~,~,info] = runAllDecoding(pPerm,verbose);
        toc
        for condPairInd = 1:length(pPerm.svm.condPairList)
            for respFeatInd = 1:length(pPerm.svm.respFeatList)
                resPerm.perfMetric(:,condPairInd,respFeatInd,permInd) = resBS{condPairInd,respFeatInd}.auc;
            end
        end
    end
end





if p.termOption.save
    disp(datestr(now))
    diary off
end
return

groupAna(p,figOption,verbose)







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
