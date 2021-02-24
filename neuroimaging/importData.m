function importData(verbose)
actuallyRun = 1;
if ~actuallyRun
    disp(['skipping ' mfilename])
    return
end
if ~exist('verbose','var')
    verbose = 1;
end

%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
subjStimList = {'jp' 'sk' 'sp' 'bm' 'sb' 'bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
    funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
        inDir  = 'a';
        outDir = 'b';
    anatPath = fullfile(repoPath,'C-derived\DecodingHR\anat\z');
    stimPath = fullfile(repoPath,'B-clean\DecodingHR\stim\160118_cyclicStim\data');
%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'funPath' 'anatPath' 'stimPath' 'inDir' 'outDir'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end
maskLabel = 'v1'; % 'v1v2v3'



for subjInd = 1:length(subjList)
    %% Get data and design
    clearvars -except tstart mask subjInd smLevel subjStimList subjList maskLabel matFun repo funPath anatPath stimPath inDir outDir noMovement runInd figOption verbose
    subj = subjList{subjInd}; subjStim = subjStimList{subjInd};
    
    switch getenv('OS')
        case 'Linux'
            funData_folderIN = fullfile('/mnt/hgfs/work/projects/160707_HRdecodingGLMdenoise/C_processing_smooth/',subj,'/run1a_preprocessing');
            labelPath = '/mnt/hgfs/work/projects/160707_HRdecodingGLMdenoise/B_acquisition/160118_cyclicStim/data';
        otherwise
            %                     data_folder = fullfile('C:\Users\Sebastien\OneDrive - McGill University\work\projects\170210_HRdecoding\C_processing\',subj,'\run1a_preprocessing');
            %                     labelDir = 'C:\Users\Sebastien\OneDrive - McGill University\work\projects\170210_HRdecoding\B_acquisition\160118_cyclicStim\data';
            funData_folderIN = fullfile(funPath,inDir,subj);
            funData_folderOUT = fullfile(funPath,outDir);
            anatData_folder = fullfile(anatPath,subj);
            labelPath = fullfile(stimPath,subjStim);
    end
    if verbose
        disp([subj ': importing from ''' funData_folderIN ''''])
    else
        disp([subj ': importing'])
    end
    
    tmp = dir(fullfile(funData_folderIN,['trun*_preprocessed.nii.gz']));
    if isempty(tmp)
        tmp = dir(fullfile(funData_folderIN,['trun*_preprocessed.nii']));
    end
    for i = 1:length(tmp)
        files{i,:} = tmp(i).name;
    end
    
    
    %Get run labels
    runList = dir(fullfile(labelPath,'*.mat'));
    label = [];
    i=1;
    while i <= length(runList)
        tmp = strsplit(runList(i).name,'__');
        if length(tmp)==3
            label = [label str2num(tmp{3}(18:end-4))];
        end
        i = i+1;
    end
    
    
    %% Load mask
    if strcmp(maskLabel,'v1v2v3')
        mask = load_nii(fullfile(anatData_folder,'v1.nii.gz'));
        tmpMask = mask.img;
        mask = load_nii(fullfile(anatData_folder,'v2.nii.gz'));
        tmpMask(logical(mask.img)) = 1;
        mask = load_nii(fullfile(anatData_folder,'v3.nii.gz'));
        tmpMask(logical(mask.img)) = 1;
        mask = tmpMask; clear tmpMask
        mask = double(flipdim(permute(mask,[3 1 2 4]),1));
    else
        maskFile = dir(fullfile(anatData_folder,[maskLabel '.nii.gz']));
        if isempty(maskFile)
            maskFile = dir(fullfile(anatData_folder,[maskLabel '.nii']));
        end
        maskFile = fullfile(anatData_folder,maskFile.name);
        mask = load_nii(maskFile);
        mask = double(flipdim(permute(mask.img,[3 1 2 4]),1));
    end
    % Censor first and last slices
    mask(:,:,[1 end]) = false;
    any3 = repmat(any(mask,3),[           1            1 size(mask,3)]);
    any2 = repmat(any(any3,2),[           1 size(any3,2)            1]);
    any1 = repmat(any(any3,1),[size(any3,1)            1            1]);
    cropMask = any1&any2; clear any1 any2 any3
    cropMask(:,:,[1 end]) = false;
    
    %% Process for all conditions
    sessionLabel = cell(3,length(files)/3);
    data = cell(3,length(files)/3);
    design = cell(3,length(files)/3);
    extraRegr = cell(3,length(files)/3);
    runInd = nan(3,length(files)/3);
    labelList = [45 135 999];
    condCount = zeros(1,3);
    for cond = 1:3
        clearvars -except tstart mask cropMask subjInd subjList subjStimList smLevel subj funData_folderIN funData_folderOUT labelDir files label labelList cond sessionLabel data design extraRegr labelList curLabel condCount sessionLabel maskLabel matFun repo funPath anatPath stimPath inDir outDir noMovement runInd figOption verbose
        close all
        curLabel = labelList(cond);
        
        %Extract runs
        for i = 1:length(files) %[1:6 19:24]
            if label(i)==curLabel
                condCount(curLabel==labelList) = condCount(curLabel==labelList)+1;
                curInd1 = cond;
                curInd2 = condCount(curLabel==labelList);
                runInd(curInd1,curInd2) = i;
                if verbose>1
                    display(['Loading condition ' num2str(cond) '; run ' num2str(curInd2) '; ' subjList{subjInd}])
                end
                %Load data
                nDataPtsToCensor = 0;
                sessionLabel{curInd1,curInd2} = str2num(files{i}(5));
                curRun = load_nii(fullfile(funData_folderIN,files{i}));
                allTr(i) = curRun.original.hdr.dime.pixdim(5);
                data{curInd1,curInd2} = double(flipdim(permute(curRun.img(:,:,:,(1+nDataPtsToCensor):end),[3 1 2 4]),1));
                %Crop data
                data{curInd1,curInd2} = data{curInd1,curInd2}(:,:,squeeze(any(any(cropMask,1),2)),:);
                data{curInd1,curInd2} = data{curInd1,curInd2}(:,squeeze(any(any(cropMask,1),3)),:,:);
                data{curInd1,curInd2} = data{curInd1,curInd2}(squeeze(any(any(cropMask,2),3)),:,:,:);
                
                %                         data{curInd1,curInd2} = data{curInd1,curInd2}(:,:,11,:);
                %                         data{curInd1,curInd2} = data{curInd1,curInd2}(1:2,1:2,1:2,:);
                %                         data{curInd1,curInd2} = data{curInd1,curInd2}(1:50,48:49,10:11,:); data{curInd1,curInd2}(1:45,:,:,:) = nan;
                %                     data{end+1} = flipdim(permute(curRun.img,[3 1 2 4]),1);
                %Gen design
                design{curInd1,curInd2} = zeros(size(data{curInd1,curInd2},4),1);
                design{curInd1,curInd2}((0:12:(120-nDataPtsToCensor)-12)+1) = 1;
                %Load extra regressors
                extraRegr{curInd1,curInd2} = dlmread(fullfile(funData_folderIN,[files{i}(1:7) '_mcParam_all.1D']));
                %                         extraRegr{end} = extraRegr{end}(1+nDataPtsToCensor):end,1:6); % Keep only movement
                extraRegr{curInd1,curInd2} = extraRegr{curInd1,curInd2}((1+nDataPtsToCensor):end,[1:6 13:18]); % Keep only Movement and derivatives; exclude movement squared
                extraRegr{curInd1,curInd2} = extraRegr{curInd1,curInd2} - repmat(mean(extraRegr{curInd1,curInd2},1),size(extraRegr{curInd1,curInd2},1),1); % Subtract the mean (will not change anything, but make the design matrix reflect better what is actually regressed out)
                for ii = 1:size(extraRegr{curInd1,curInd2},2)
                    extraRegr{curInd1,curInd2}(:,ii) = extraRegr{curInd1,curInd2}(:,ii)./max(abs(extraRegr{curInd1,curInd2}(:,ii))); % Scale between zero and one
                end
            end
        end
    end
    %             tstart = tic;
    
    stimdur = 6;
    tr = round(curRun.original.hdr.dime.pixdim(5)*10)/10;
    
    
    % Check various constants
    if verbose>1
        fprintf('There are %d runs in total.\n',length(design));
        fprintf('The dimensions of the data for the first run are %s.\n',mat2str(size(data{1})));
        fprintf('The stimulus duration is %.6f seconds.\n',stimdur);
        fprintf('The sampling rate (TR) is %.6f seconds.\n',tr);
    end
    
    clearvars -except tstart mask cropMask subjInd smLevel subjStimList subjList subj funData_folderIN funData_folderOUT labelDir data design extraRegr sessionLabel sessModel stimdur tr maskLabel repo funPath anatPath stimPath inDir outDir noMovement runInd figOption verbose
    
    
    %% Rorganize data
    d.runInd = runInd';
    d.repLabel = (1:size(d.runInd,1))'*ones(1,size(d.runInd,2));
    d.condLabel = ones(size(d.runInd,1),1)*[1 2 3];
    d.sessLabel = cell2mat(sessionLabel)';
    d.data = data'; clear data
    d.design = design'; clear design
    d.extraRegr = extraRegr'; clear extraRegr
    
    [~,b] = sort(d.runInd(:));
    fieldList = fields(d);
    for fieldInd = 1:length(fieldList)
        d.(fieldList{fieldInd}) = d.(fieldList{fieldInd})(b);
    end
    sessLabel = d.sessLabel;
    for fieldInd = 1:length(fieldList)
        d1.(fieldList{fieldInd}) = d.(fieldList{fieldInd})(sessLabel==1);
        d2.(fieldList{fieldInd}) = d.(fieldList{fieldInd})(sessLabel==2);
        d.(fieldList{fieldInd}) = [];
    end
    clear d
    d2.repLabel = d2.repLabel - min(d2.repLabel) + 1;
    
    d.fun = rmfield(cat(2,d1,d2),{'sessLabel' 'runInd'}); clear d1 d2
    d.info = '1 X sess';
    p.tr = 1;
    p.stimDur = stimdur;
    p.doMotion = 0;
    p.masks.roiMasks.(maskLabel) = logical(mask);
    p.masks.cropMask = logical(cropMask);
    
    %% Save
    disp([subj ': saving croped data'])
    if ~exist(funData_folderOUT,'dir')
        mkdir(funData_folderOUT)
    end
    outName = fullfile(funData_folderOUT,subj);
    save([outName '.mat'],'d','p','-v7.3');
    disp([subj ': saved to ''' outName ''''])
end