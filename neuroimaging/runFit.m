function runFit

actuallyRun = 0;
noMovement = 1;
saveFig = 0;
exclude = 1;
exclusion.subj = {'03sk'};
exclusion.sess = 1;
exclusion.run = 4;
exclusion.cond = 2;


if ismac
    repo = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repo = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
    funDir = 'C-derived\DecodingHR\fun';
        inDir = 'x';
        outDir = 'y';
    anatDir = 'C-derived\DecodingHR\anat\z';
    stimDir = 'B-clean\DecodingHR\stim\160118_cyclicStim\data';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repo' 'funDir' 'anatDir' 'stimDir'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end

% maskLabel = 'v1v2v3';
maskLabel = 'v1';

% subjList = {'03sk'};
% subjStimList = {'sk'};
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
subjStimList = {'jp' 'sk' 'sp' 'bm' 'sb' 'bj'};


disp('Please get the data from the repo if not done already')
disp(['data repo: ' repo])
disp(['IN: anatomical V1 roi (' fullfile(anatDir) ')'])
disp(['IN: preprocessed functionals (' fullfile(funDir,inDir) ')'])
disp(['IN: stimulus timing (' fullfile(stimDir) ')'])
disp('F(IN)=OUT: 2-df sinusoidal fit to single voxel time series')
disp(['OUT: fit params and stats + HRF estimates (' fullfile(funDir,outDir) ')'])

if ~actuallyRun
    disp('Not actually running because way too long')
    
    subjInd = 1;
    
    % Random Effect Model
    f = figure('WindowStyle','docked');
    
    subplot(1,16+25,25+(1:16))
    if noMovement
        load(fullfile(repo,funDir,outDir,subjList{subjInd},'v1SinCos_1perRun.mat'),'results')
        run1 = results.OLS.mixed.designmatrix(:,[1 2 67:68]);
    else
        load(fullfile(repo,funDir,outDir,subjList{subjInd},'v1SinCos_1perRun_move12.mat'),'results')
        run1 = results.OLS.mixed.designmatrix(:,[1 2 67:80]);
    end
    %scale for max contrast matix
    tmp = run1(any(run1,2),:);
    tmp(:,3) = max(tmp(:,1));
    tmp(:,4) = tmp(:,4)./max(abs(tmp(:,4)));
    imagesc(tmp); colormap gray
    title('Sinusoidal Response');
    ax1 = gca; ax1.XTick = []; ax1.YTick = [];
%     xlabel('Regressors');
    ylabel('TRs');
    ax1.XTick = 1:size(run1,2);
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    ax1.XTickLabel = {'sin' 'cos' 'constant' 'drift' 'x' 'y' 'z' 'pitch' 'roll' 'yaw' 'x''' 'y''' 'z''' 'pitch''' 'roll''' 'yaw'''};
    ax1.XTickLabelRotation = -90;
    ax1.YAxis.Color = 'none';
    ax1.YAxis.Label.Visible = 'on';
    ax1.YAxis.Label.Color = 'k';
    
    xVal = diff(ax1.XLim);
    xWidth = ax1.Position(3);
    
    subplot(1,16+25,0+(1:25))
    if noMovement
        tmp = load(fullfile(repo,funDir,outDir,subjList{1},'v1resp_1perRun_resp.mat'));
        run1 = tmp.results.OLS.mixed.designmatrix(:,[1:12 397]);
    else
        tmp = load(fullfile(repo,funDir,outDir,subjList{1},'v1resp_1perRun_move12_resp.mat'));
        run1 = tmp.results.OLS.mixed.designmatrix(:,[1:12 397:409]);
    end
    tmp = run1(any(run1,2),:);
    tmp(:,end) = tmp(:,end)./max(abs(tmp(:,end)));
    imagesc(tmp); colormap gray
    title('Model-Free Response');
    ax2 = gca; ax2.XTick = []; ax2.YTick = [];
%     xlabel('Regressors');% ylabel('TRs');
    ax2.XTick = 1:size(run1,2);
    ax2.Box = 'off';
    ax2.TickDir = 'out';
    ax2.XTickLabel = cat(1,cellstr([repmat('t+',12,1) num2str((0:11)','%-d')]),{'drift'},{'x' 'y' 'z' 'pitch' 'roll' 'yaw' 'x''' 'y''' 'z''' 'pitch''' 'roll''' 'yaw'''}');
    ax2.XTickLabelRotation = -90;
    ax2.YAxis.Color = 'none';
    ax2.YAxis.Label.Visible = 'on';
    ax2.YAxis.Label.Color = 'k';
    
    suptitle('Design Matrices (random-effect)');
    
    ax1.XAxis.FontSize = ax1.XAxis.FontSize*0.8;
    ax2.XAxis.FontSize = ax2.XAxis.FontSize*0.8;
    
    ax = [ax1 ax2];
    for i = 1:length(ax)
        curAxe = ax(i);
        curAxeB = axes();
        linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
        curAxeB.Color = 'none'; curAxeB.Box = 'on';
        curAxeB.XTick = []; curAxeB.YTick = [];
%         delta = mode(diff(curAxe.XTick));
%         curAxeB.XTick = [curAxe.XTick curAxe.XTick(end)+delta]-delta/2;
%         curAxeB.TickLength = [0 0];
%         curAxeB.XTickLabel = [];
%         curAxeB.XGrid = 'on';
%         curAxeB.GridAlpha = 1;
    end
    if saveFig
        saveas(f,fullfile(repo,funDir,outDir,'designMatrices_randomEffect'))
        disp(fullfile(repo,funDir,outDir,'designMatrices_randomEffect.fig'))
    end
    
    
    % Fixed Effect model
    f = figure('WindowStyle','docked');
    if noMovement
        load(fullfile(repo,funDir,outDir,subjList{subjInd},'v1SinCos_1perRun.mat'),'results')
%         run1 = results.OLS.fixed.designmatrix(:,[1 2 67:68]);
    else
        load(fullfile(repo,funDir,outDir,subjList{subjInd},'v1SinCos_1perRun_move12.mat'),'results')
%         run1 = results.OLS.fixed.designmatrix(:,[1 2 67:80]);
    end
    sessInd = results.OLS.fixed.designmatrixSessInd==1;
    run1 = results.OLS.fixed.designmatrix(sessInd,:);
    runOrderLabel = results.OLS.fixed.designmatrixRunInd(sessInd,:);
    run1 = run1(:,~all(run1==0,1));
    runOrderLabelList = unique(runOrderLabel);
    tmp = [];
    for i = 1:length(runOrderLabelList)
        trInd = runOrderLabel==runOrderLabelList(i);
        tmp = cat(1,tmp,run1(trInd,:));
    end
    tmp(:,7:end) = run1(:,7:end);
    run1 = tmp;
    tmp = run1(:,7:2:end);
    run1(:,7:2:end) = run1(:,7:2:end)./mode(tmp(tmp~=0));
    tmp = run1(:,8:2:end);
    run1(:,8:2:end) = run1(:,8:2:end)./max(abs(tmp(tmp~=0)));
    
    imagesc(run1); colormap gray
    title({'Fixed-Effect Model' '\fontsize{8}Sinusoidal Response'});
    ax1 = gca; ax1.XTick = []; ax1.YTick = [];
%     xlabel('Regressors');
    ylabel('TRs');
    ax1.XTick = [1.5 3.5 5.5];
    ax1.Box = 'off';
    ax1.TickDir = 'out';
    ax1.XTickLabel = {'ori1' 'ori2' 'plaid'};
    ax1.XTickLabelRotation = -90;
    ax1.YAxis.Color = 'none';
    ax1.YAxis.Label.Visible = 'on';
    ax1.YAxis.Label.Color = 'k';
    
    ax = [ax1];
    for i = 1:length(ax)
        curAxe = ax(i);
        curAxeB = axes();
        linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
        curAxeB.Color = 'none'; curAxeB.Box = 'on';
        curAxeB.XTick = []; curAxeB.YTick = [];
%         delta = mode(diff(curAxe.XTick));
%         curAxeB.XTick = [curAxe.XTick curAxe.XTick(end)+delta]-delta/2;
%         curAxeB.TickLength = [0 0];
%         curAxeB.XTickLabel = [];
%         curAxeB.XGrid = 'on';
%         curAxeB.GridAlpha = 1;
    end
    curAxeB.XTick = mean([6.5 size(run1,2)+0.5]);
    curAxeB.XTickLabel = {'constant+drift for each run'};
    curAxeB.TickDir = 'out';
    
    if saveFig
        saveas(f,fullfile(repo,funDir,outDir,'designMatrices_fixedEffect'))
        disp(fullfile(repo,funDir,outDir,'designMatrices_fixedEffect.fig'))
    end
    return
end

    

for smLevel = {''}
% for smLevel = {'' '_sm1.50' '_sm2.00' '_sm3.00' '_sm4.00' '_sm1.00' '_sm1.25' '_sm1.75' '_sm2.50'}
    for subjInd = 1:length(subjList)
%         try
            
            %% Get data and design
            clearvars -except tstart mask subjInd smLevel subjStimList subjList maskLabel matFun repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion
            subj = subjList{subjInd}; subjStim = subjStimList{subjInd};
            
            switch getenv('OS')
                case 'Linux'
                    funData_folderIN = fullfile('/mnt/hgfs/work/projects/160707_HRdecodingGLMdenoise/C_processing_smooth/',subj,'/run1a_preprocessing');
                    labelDir = '/mnt/hgfs/work/projects/160707_HRdecodingGLMdenoise/B_acquisition/160118_cyclicStim/data';
                otherwise
%                     data_folder = fullfile('C:\Users\Sebastien\OneDrive - McGill University\work\projects\170210_HRdecoding\C_processing\',subj,'\run1a_preprocessing');
%                     labelDir = 'C:\Users\Sebastien\OneDrive - McGill University\work\projects\170210_HRdecoding\B_acquisition\160118_cyclicStim\data';
                    funData_folderIN = fullfile(repo,funDir,inDir,subj);
                    funData_folderOUT = fullfile(repo,funDir,outDir,subj);
                    if ~exist(funData_folderOUT,'dir'); mkdir(funData_folderOUT); end
                    anatData_folder = fullfile(repo,anatDir,subj);
                    labelDir = fullfile(repo,stimDir);
            end
            
            tmp = dir(fullfile(funData_folderIN,['trun*_preprocessed' smLevel{1} '.nii.gz']));
            if isempty(tmp)
                tmp = dir(fullfile(funData_folderIN,['trun*_preprocessed' smLevel{1} '.nii']));
            end
            for i = 1:length(tmp)
                files{i,:} = tmp(i).name;
            end
            
            
            
            %Get run labels
            runList = dir(fullfile(labelDir,subjStim,'*.mat'));
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
            any3 = repmat(any(mask,3),[           1            1 size(mask,3)]);
            any2 = repmat(any(any3,2),[           1 size(any3,2)            1]);
            any1 = repmat(any(any3,1),[size(any3,1)            1            1]);
            mask = any1&any2; clear any1 any2 any3
            
            %% Process for all conditions
            sessionLabel = cell(3,length(files)/3);
            data = cell(3,length(files)/3);
            design = cell(3,length(files)/3);
            extraRegr = cell(3,length(files)/3);
            runInd = nan(3,length(files)/3);
            labelList = [45 135 999];
            condCount = zeros(1,3);
            for cond = 1:3
                clearvars -except tstart mask subjInd subjList subjStimList smLevel subj funData_folderIN funData_folderOUT labelDir files label labelList cond sessionLabel data design extraRegr labelList curLabel condCount sessionLabel maskLabel matFun repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion
                close all
                curLabel = labelList(cond);
                
                %Extract runs
                for i = 1:length(files) %[1:6 19:24]
                    if label(i)==curLabel
                        condCount(curLabel==labelList) = condCount(curLabel==labelList)+1;
                        curInd1 = cond;
                        curInd2 = condCount(curLabel==labelList);
                        runInd(curInd1,curInd2) = i;
                        display(['Loading condition ' num2str(cond) '; run ' num2str(curInd2) '; ' subjList{subjInd}])
                        %Load data
                        nDataPtsToCensor = 12;
                        sessionLabel{curInd1,curInd2} = str2num(files{i}(5));
                        curRun = load_nii(fullfile(funData_folderIN,files{i}));
                        allTr(i) = curRun.original.hdr.dime.pixdim(5);
                        data{curInd1,curInd2} = double(flipdim(permute(curRun.img(:,:,:,(1+nDataPtsToCensor):end),[3 1 2 4]),1));
                        %Crop data
                        data{curInd1,curInd2} = data{curInd1,curInd2}(:,:,squeeze(any(any(mask,1),2)),:);
                        data{curInd1,curInd2} = data{curInd1,curInd2}(:,squeeze(any(any(mask,1),3)),:,:);
                        data{curInd1,curInd2} = data{curInd1,curInd2}(squeeze(any(any(mask,2),3)),:,:,:);

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
            fprintf('There are %d runs in total.\n',length(design));
            fprintf('The dimensions of the data for the first run are %s.\n',mat2str(size(data{1})));
            fprintf('The stimulus duration is %.6f seconds.\n',stimdur);
            fprintf('The sampling rate (TR) is %.6f seconds.\n',tr);
            
            
            clearvars -except tstart mask subjInd smLevel subjStimList subjList subj funData_folderIN funData_folderOUT labelDir data design extraRegr sessionLabel sessModel stimdur tr maskLabel repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion
            %% GLMdenoise on all sessions (not split)
            % - - - - - - - - - -
            % o x ------------
            splitIn = 1;
            display([subj '; split in ' num2str(splitIn)])
            
            splitPossible = 8;
            splitData = data;
            splitDesign = design;
            splitExtraRegr = extraRegr;
            splitSessionLabel = sessionLabel;
            splitRunInd = num2cell(runInd);
            
            splitData = cat(2,splitData(1,:),splitData(2,:),splitData(3,:)); splitData(cellfun('isempty',splitData)) = [];
            splitDesign = cat(2,splitDesign(1,:),splitDesign(2,:),splitDesign(3,:)); splitDesign(cellfun('isempty',splitDesign)) = [];
            splitExtraRegr = cat(2,splitExtraRegr(1,:),splitExtraRegr(2,:),splitExtraRegr(3,:)); splitExtraRegr(cellfun('isempty',splitExtraRegr)) = [];
            splitSessionLabel = cat(2,splitSessionLabel(1,:),splitSessionLabel(2,:),splitSessionLabel(3,:)); splitSessionLabel(cellfun('isempty',splitSessionLabel)) = [];
            splitRunInd = cat(2,splitRunInd(1,:),splitRunInd(2,:),splitRunInd(3,:)); splitRunInd(cellfun('isempty',splitRunInd)) = [];
            
            splitRun_tmp = cell(length(splitData),splitIn);
            splitInd_tmp = cell(length(splitData),splitIn);
            splitData_tmp = cell(length(splitData),splitIn);
            splitDesign_tmp = cell(length(splitData),splitIn);
            splitExtraRegr_tmp = cell(length(splitData),splitIn);
            splitSessionLabel_tmp = cell(length(splitData),splitIn);
            splitRunInd_tmp = cell(length(splitData),splitIn);
            for i = 1:length(splitData)
                %remove another cycle at the begining, keep the end
                splitData{i}(:,:,:,1:12) = [];
                splitDesign{i}(1:12) = [];
                splitExtraRegr{i}(1:12,:) = [];
                
                for ii = 1:splitIn
                    splitRun_tmp{i,ii} = i;
                    splitInd_tmp{i,ii} = 1+(ii-1)*(splitPossible/splitIn)*12:(ii)*(splitPossible/splitIn)*12;
                    splitData_tmp{i,ii} = splitData{i}(:,:,:,splitInd_tmp{i,ii});
                    splitDesign_tmp{i,ii} = splitDesign{i}(splitInd_tmp{i,ii});
                    splitExtraRegr_tmp{i,ii} = splitExtraRegr{i}(splitInd_tmp{i,ii},:);
                    splitSessionLabel_tmp{i,ii} = splitSessionLabel{i};
                    splitRunInd_tmp{i,ii} = splitRunInd{i};
                end
            end
            splitRun = cell(1,numel(splitRun_tmp));
            splitInd = cell(1,numel(splitInd_tmp));
            splitData = cell(1,numel(splitData_tmp));
            splitDesign = cell(1,numel(splitDesign_tmp));
            splitExtraRegr = cell(1,numel(splitExtraRegr_tmp));
            splitSessionLabel = cell(1,numel(splitSessionLabel_tmp));
            splitRunInd = cell(1,numel(splitRunInd_tmp));
            for ii = 1:splitIn
                splitRun(ii:splitIn:numel(splitRun_tmp)) = splitRun_tmp(:,ii);
                splitInd(ii:splitIn:numel(splitInd_tmp)) = splitInd_tmp(:,ii);
                splitData(ii:splitIn:numel(splitData_tmp)) = splitData_tmp(:,ii);
                splitDesign(ii:splitIn:numel(splitDesign_tmp)) = splitDesign_tmp(:,ii);
                splitExtraRegr(ii:splitIn:numel(splitExtraRegr_tmp)) = splitExtraRegr_tmp(:,ii);
                splitSessionLabel(ii:splitIn:length(splitDesign)) = splitSessionLabel_tmp(:,ii);
                splitRunInd(ii:splitIn:length(splitDesign)) = splitRunInd_tmp(:,ii);
            end
            clearvars -except tstart mask subjInd smLevel subjStimList subjList subj funData_folderIN funData_folderOUT labelDir data design extraRegr sessionLabel sessModel splitDesign splitData splitIn splitExtraRegr splitSessionLabel stimdur tr maskLabel repo funDir anatDir stimDir inDir outDir noMovement runInd splitRunInd exclude exclusion
            
            
            if noMovement
                disp('--- Response Time Course Extraction ---')
                ana = 'resp';
                [results] = GLMresp(splitDesign,splitData,stimdur,tr,ana,[],struct('sessionLabel',{splitSessionLabel},'splitedIn',splitIn),splitIn);
                outName = fullfile(funData_folderOUT,[maskLabel ana '_' num2str(splitIn) 'perRun' smLevel{1}]);
                save([outName '_resp.mat'],'results'); clear results
                
                disp('--- Response Sinusoidal Fit ---')
                ana = 'SinCos';
                outName = fullfile(funData_folderOUT,[maskLabel ana '_' num2str(splitIn) 'perRun' smLevel{1}]);
                [resultsTmp,~] = GLMsinCos3(splitDesign,splitData,stimdur,tr,ana,[],struct('sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn);
                resultsTmp.OLS.fixed = []; resultsTmp.OLS.fixed_sessReg = []; resultsTmp.OLS.mixed_sessRm = [];
                if exclude && ismember(subjList(subjInd),exclusion.subj)
                    exclInd = ismember(exclusion.subj,subjList(subjInd));
                    curExclusion.sess = exclusion.sess(exclInd);
                    curExclusion.run = exclusion.run(exclInd);
                    curExclusion.cond = exclusion.cond(exclInd);
                    [results,~] = GLMsinCos5(splitDesign,splitData,stimdur,tr,ana,[],struct('sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn,curExclusion);
                else
                    [results,~] = GLMsinCos5(splitDesign,splitData,stimdur,tr,ana,[],struct('sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn);
                end
                results.OLS.mixed = resultsTmp.OLS.mixed; clear resultsTmp
                results.mask = mask;
                save([outName '.mat'],'results','-v7.3');
            else
                ana = 'resp';
                [results] = GLMresp(splitDesign,splitData,stimdur,tr,ana,[],struct('extraregressors',{splitExtraRegr},'sessionLabel',{splitSessionLabel},'splitedIn',splitIn),splitIn);
                outName = fullfile(funData_folderOUT,[maskLabel ana '_' num2str(splitIn) 'perRun_move12' smLevel{1}]);
                save([outName '_resp.mat'],'results'); clear results
                
                ana = 'SinCos';
                outName = fullfile(funData_folderOUT,[maskLabel ana '_' num2str(splitIn) 'perRun_move12' smLevel{1}]);
                [resultsTmp,~] = GLMsinCos3(splitDesign,splitData,stimdur,tr,ana,[],struct('extraregressors',{splitExtraRegr},'sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn);
                resultsTmp.OLS.fixed = []; resultsTmp.OLS.fixed_sessReg = []; resultsTmp.OLS.mixed_sessRm = [];
                if exclude && ismember(subjList(subjInd),exclusion.subj)
                    exclInd = ismember(exclusion.subj,subjList(subjInd));
                    curExclusion.sess = exclusion.sess(exclInd);
                    curExclusion.run = exclusion.run(exclInd);
                    curExclusion.cond = exclusion.cond(exclInd);
                    [results,~] = GLMsinCos5(splitDesign,splitData,stimdur,tr,ana,[],struct('extraregressors',{splitExtraRegr},'sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn,curExclusion);
                else
                    [results,~] = GLMsinCos5(splitDesign,splitData,stimdur,tr,ana,[],struct('extraregressors',{splitExtraRegr},'sessionLabel',{splitSessionLabel},'runLabel',{splitRunInd},'splitedIn',splitIn),splitIn);
                end
                results.OLS.mixed = resultsTmp.OLS.mixed; clear resultsTmp
                results.mask = mask;
                results.inputs.opt.runLabel = splitRunInd;
                save([outName '.mat'],'results','-v7.3');
            end

%             [results,dataDetrend] = GLMsinCos5(splitDesign,splitData,stimdur,tr,ana,[],struct('extraregressors',{splitExtraRegr},'sessionLabel',{splitSessionLabel},'splitedIn',splitIn),splitIn);
%             save([outName '_detrend.mat'],'dataDetrend','-v7.3');
%             clear results dataDetrend
            
            disp([subj '; split in ' num2str(splitIn) ': done'])
            
            telapsed = toc;
            disp(['Whole thing took: ' num2str(telapsed) 'sec'])
    end
end
