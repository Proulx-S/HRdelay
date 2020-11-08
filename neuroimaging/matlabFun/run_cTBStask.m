function [ output_args ] = run_cTBStask(allFile)

for taskInd = 1:length(allFile)
    clearvars -except taskInd allFile
    display(['Performing ' char(allFile{taskInd})])
    load(allFile{taskInd})
    
    
    
    [a,b,c] = fileparts(dataFileIn);
    load(fullfile(dataDirIn,[b '.mat']),'metric')
    
    tmp.runList = metric.runList;
    tmp.run = metric.run;
    tmp.trial = metric.(['trial' num2str(p.trialsPerObs)]);
    metric = tmp; clear tmp
    
    %Compile it
    runList = dir([labelDir '/*.mat']);
    label = [];
    i=1;
    while i <= length(runList)
        tmp = strsplit(runList(i).name,'__');
        if length(tmp)==3
            label = [label str2num(tmp{3}(18:end-4))];
        end
        i = i+1;
    end
    %         switch dataType
    %             case 'trial'
    %                 metricCompiled = compileRuns(metric,label,fitType,dataType);
    %                 metricCompiled_thresh = compileRuns(metric,label,fitType,'run');
    %             case 'run'
    metricCompiled = compileRuns(metric,label,fitType,dataType);
    metricCompiled_thresh = metricCompiled;
    %         end
    %Threshold it
    %based on r^2
    
    if p.thresholdData
        switch p.thresholdP
            case 0.05
                threshRsquared = 0.0358;
            case 0.01
                threshRsquared = 0.0610;
                caG2G1se 0.001
                threshRsquared = 0.0975;
            case 0.0001
                threshRsquared = 0.1336;
            otherwise
                warning('threshRsquared not precompiled for requested p, will estimate it. Should be fine as long as you have very large noisy data set')
                threshRsquared = findRsqaredAtP(metricCompiled,p.thresholdP);
        end
        metricCompiled = applyThresh(metricCompiled,metricCompiled_thresh,'Rsquared',threshRsquared);%0.2435^2
    end
    metric = metricCompiled; clear metricCompiled
    %based on phase locking
    
    if p.PL
        metric = applyPhaseLockingThresh(metric);
    end
    
    
    %% Extract baseline data
    %sanity check
    if find(diff(metric.sessionLabel,[],3))
        error
    end
    
    baseInd = metric.sessionLabel(1,:,1)==1;
    
    %amp
    baseData = metric.ampRect(:,baseInd,:);
    baseData_mean = mean(baseData,2);
    baseData_std = std(baseData,[],2);
    
    baseDataG = (metric.ampRect(:,baseInd,2) + metric.ampRect(:,baseInd,1))/2;
    baseDataG_mean = mean(baseDataG,2);
    baseDataG_std = std(baseDataG,[],2);
    
    baseDataG2G1 = metric.ampRect(:,baseInd,2) - metric.ampRect(:,baseInd,1);
    baseDataG2G1_mean = mean(baseDataG2G1,2);
    baseDataG2G1_std = std(baseDataG2G1,[],2);
    
    baseDataPG = metric.ampRect(:,baseInd,3) - baseDataG;
    baseDataPG_mean = mean(baseDataPG,2);
    baseDataPG_std = std(baseDataPG,[],2);
    
    baseDataPG1 = metric.ampRect(:,baseInd,3) - metric.ampRect(:,baseInd,1);
    baseDataPG1_mean = mean(baseDataPG1,2);
    baseDataPG1_std = std(baseDataPG1,[],2);
    
    baseDataPG2 = metric.ampRect(:,baseInd,3) - metric.ampRect(:,baseInd,2);
    baseDataPG2_mean = mean(baseDataPG2,2);
    baseDataPG2_std = std(baseDataPG2,[],2);
    
    
    
    %% Extract post data
    postInd = metric.sessionLabel(1,:,1)==2;
    
    %amp
    postData = metric.ampRect(:,postInd,:);
    postData_z = postData-repmat(baseData_mean,[1 length(find(postInd)) 1]);
    postData_z = postData_z./repmat(baseData_std,[1 length(find(postInd)) 1]);
    
    postDataG1_z = postData_z(:,:,1);
    postDataG2_z = postData_z(:,:,2);
    postDataP_z = postData_z(:,:,3);
    
    postDataG = (metric.ampRect(:,postInd,2) + metric.ampRect(:,postInd,1))/2;
    postDataG_z = postDataG-repmat(baseDataG_mean,[1 length(find(postInd)) 1]);
    postDataG_z = postDataG_z./repmat(baseDataG_std,[1 length(find(postInd)) 1]);
    
    postDataG2G1 = metric.ampRect(:,postInd,2) - metric.ampRect(:,postInd,1);
    postDataG2G1_z = postDataG2G1-repmat(baseDataG2G1_mean,[1 length(find(postInd)) 1]);
    postDataG2G1_z = postDataG2G1_z./repmat(baseDataG2G1_std,[1 length(find(postInd)) 1]);
    
    postDataPG = metric.ampRect(:,postInd,3) - postDataG;
    postDataPG_z = postDataPG-repmat(baseDataPG_mean,[1 length(find(postInd)) 1]);
    postDataPG_z = postDataPG_z./repmat(baseDataPG_std,[1 length(find(postInd)) 1]);
    
    postDataPG1 = metric.ampRect(:,postInd,3) - metric.ampRect(:,postInd,1);
    postDataPG1_z = postDataPG1-repmat(baseDataPG1_mean,[1 length(find(postInd)) 1]);
    postDataPG1_z = postDataPG1_z./repmat(baseDataPG1_std,[1 length(find(postInd)) 1]);
    
    postDataPG2 = metric.ampRect(:,postInd,3) - metric.ampRect(:,postInd,2);
    postDataPG2_z = postDataPG2-repmat(baseDataPG2_mean,[1 length(find(postInd)) 1]);
    postDataPG2_z = postDataPG2_z./repmat(baseDataPG2_std,[1 length(find(postInd)) 1]);
    
    %% Output to nii
    
    maskName = fullfile(subj,'run1c_masking','V1.nii.gz');
    [err, V, Info, ErrMessage] = BrikLoad (maskName);
    
    dataList = {'postDataG1_z' 'postDataG2_z' 'postDataP_z' 'postDataG_z' 'postDataPG_z' 'postDataPG1_z' 'postDataPG2_z' 'postDataG2G1_z'};
    out = fullfile('images',subj);
    mkdir(out)
    
    av = 'all'; % 'all' 'non' 'half'
    for i = 1:length(dataList)
        eval(['curData = ' dataList{i} ';']);
        tPts = size(curData,2);
        switch av
            case 'all'
                curIm = nan([size(V) 1]);
            case 'none'
                curIm = nan([size(V) tPts]);
            case 'half'
                curIm = nan([size(V) 2]);
        end
        
        for vox = 1:size(metric.voxInd,1)
            switch av
                case 'all'
                    tmp = mean(curData(vox,:),2);
                case 'none'
                    tmp = curData(vox,:);
                case 'half'
                    tmp = cat(2,mean(curData(vox,1:ceil(tPts/2)),2),mean(curData(vox,ceil(tPts/2)+1:end),2));
            end
            curIm(metric.voxInd(vox,1),metric.voxInd(vox,2),metric.voxInd(vox,3),:) = tmp;
        end
        writeNIFTI(fullfile(out,dataList{i}(9:end)), curIm, Info)
    end
end
