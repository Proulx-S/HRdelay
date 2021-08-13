%% Info
% This script reproduces all figures from S. Proulx et al., Neuroimage, in
% preperation.
clear all
close all

%% Initiation
anaInit; % where to define paths and some general configurations

%% Get data (preprocessed BOLD volumes, only V1 ROI voxels, one .mat file per subject)
downloadData;

%% Extract responses from timeseries
extractResponses(p);

%% Define the retinotopic representation of the stimulus field of view (data-driven with priors from a probabilistic retinotopic atlas)
doWhat = 'download';
% 'run' -> run it yourself (several minutes)
% 'download' -> download precomputed data from repository
% 'run_forced' -> same as above, but forces to rerun instead of loading
% locally saved data
% 'download_forced' -> same as above, but forces to redownload instead of loading
% locally saved data
processFov(p,doWhat);

%% Feature selection
processFeatSel(p)









%% Feature selection parameters
% Within visual field region of the stimulus
p.featSel.fov.doIt = 1;
p.featSel.fov.threshMethod = 'empirical'; % 'empirical' 'ecc'
p.featSel.fov.areaLabel = 'v1';
p.featSel.fov.threshVal = [0.75 7]; % threshMethod='ecc'
p.featSel.fov.percentile = 20; % threshMethod='ecc'
% Activated voxels (fixed-effect)
p.featSel.act.doIt = 0;
p.featSel.act.threshMethod = '%ile'; % '%ile' 'p' 'fdr'
p.featSel.act.threshVal = 0.05; % threshMethod='p' or 'fdr'
p.featSel.act.percentile = 20; % threshMethod='%ile'
% Activated voxels (random-effect)
p.featSel.respVecSig.doIt = 1;
p.featSel.respVecSig.threshMethod = 'p'; % '%ile' 'p' 'fdr'
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
p.featSel.global.method = 'allData';
% 'allData'-> featSel uses all three conditions, irrespective of the condition pairs to be decoded
% 'custom1'-> featSel of active voxels uses all three conditions but featSel of discriminant voxels uses only the conditions to be decoded
% 'custom2'-> featSel of active and most discriminant voxels uses only the conditions to be decoded

%% SVM parameters
p.svm.kernel.type = 'lin';
p.svm.complexSpace = 'bouboulisDeg1'; % 'bouboulisDeg1' 'bouboulisDeg2'
p.svm.doWithin = 0;

%% Channel parameters
p.svm.condPairList = {'grat1VSgrat2' 'grat1VSplaid' 'grat2VSplaid'};
% p.svm.respFeatList = {'cartNoDelay' 'delay' 'cart'};
p.svm.respFeatList = {'delay' 'cartNoDelay' 'cart'};


if 0
    importData(verbose)
    applyAreaMask(figOption)
    %     processWaveletResponses(figOption,verbose)
end
if 0
    extractResponses(p,figOption,verbose)
end
p.featSel.fov.resFile = processFov(p);
if 0
    processFeatSel(p)
    [resBS,resBShr,resWS,f,info,decodingOut] = runAllDecoding(p);
end
visualizeOthers(p)
if 0
    visualizeFeatSel(p)
    
end
if 0
    if ~exist('resBS','var')
        decodingOut = 'C:\Users\sebas\OneDrive - McGill University\dataBig\C-derived\DecodingHR\fun\f_2021-06-15_allData\decoding';
        load([decodingOut ''],'resBS','info');
    end
    plotAllDecoding2(p,resBS,info);
    plotAllDecoding(p,resBS,info);
    statsAllDecoding(p,resBS,info)
end
if p.perm.doIt
    if 0
        disp('***************************')
        disp('Permutation test: computing')
        disp('***************************')
        for permInd = 1:p.perm.n
            disp(['Perm ' num2str(permInd) '/' num2str(p.perm.n)])
            if permInd==1
                resAll = [];
                featSel_fov = [];
            end
            [resPall,resAll] = permuteLabels(p,resAll);
            [featSelPall,featSel_fov] = processFeatSel(p,resPall,featSel_fov);
            [resBSP,resBShrP,resWSP,fP,infoP] = runAllDecoding(p,resPall,featSelPall);
            if permInd == 1
                auc = nan([size(resBSP) p.perm.n]);
            end
            tmp = reshape([resBSP{:}],size(resBSP));
            tmp = reshape(permute([tmp.auc],[2 3 1]),[size(tmp) size(tmp(1).auc,1)]);
            auc(:,:,permInd) = mean(tmp,3); clear tmp
        end
        disp('**********************')
        disp('Permutation test: done')
        disp('**********************')
        % Save permutations
        disp('Saving permutations')
        auc = permute(auc,[3 1 2]);
        for i = 1:numel(auc(1,:))
            resBS{i}.subj.aucP = auc(:,i);
        end
        save([decodingOut 'P'],'resBS','info');
        disp(['Permutations saved to ' decodingOut 'P'])
    else
        decodingOut = 'C:\Users\sebas\OneDrive - McGill University\dataBig\C-derived\DecodingHR\fun\f_2021-06-15_allData\decoding';
        load([decodingOut 'P'],'resBS','info');
        condInd = 1:3;
        info.condPairList(condInd);
        respInd = 1:2;
        info.respFeatList(respInd);
        
        tmp = reshape([resBS{:}],size(resBS));
        tmp = reshape([tmp.subj],size(resBS));
        auc = permute(reshape([tmp.auc],[length(tmp(1).auc) size(resBS)]),[2 3 1]);
        aucP = permute(reshape([tmp.aucP],[length(tmp(1).aucP) size(resBS)]),[2 3 1]);
        
        auc = mean(auc,1);
        [H,P,CI,STATS] = ttest(squeeze(auc)',0.5,'tail','right');
        
        auc = mean(auc,3);
        aucP = mean(aucP,1);
        permP = sum(auc<aucP,3)./size(aucP,3);
        
        respInd = 1:3;
        info.respFeatList(respInd);
        disp('****Permutation test (all)****')
        for i = 1:length(respInd)
            disp(['* ' info.respFeatList{respInd(i)} ': AUC=' num2str(auc(i)) ', p=' num2str(permP(respInd(i)))])
        end
        disp('************************')
    end
    plotAllDecoding3(p,resBS,info);
end

if p.termOption.save
    disp(datestr(now))
    diary off
end
return







if 0
    chan = processChanHr(p,resBShr,info);
    f = plotChanHr(p,chan);
    statsChanHr(p,chan);
end


return



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
        t1 = toc;
        processFeatSel(pPerm)
        t2 = toc;
        [resBS,~,~,~,info] = runAllDecoding(pPerm,verbose);
        t3 = toc;
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