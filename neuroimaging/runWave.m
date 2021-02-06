function runWave(verbose,figOption)
if ~exist('verbose','var')
    verbose = 1;
end

actuallyRun = 1;
noMovement = 1;
if exist('figOption','var') && isfield(figOption,'save')
    saveFig = figOption.save;
else
    saveFig = 0;
end
exclude = 1;
exclusion.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
exclusion.subj = 2;
exclusion.sess = {1};
exclusion.run = {4};
exclusion.cond = {2};


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


if ~actuallyRun
    filename = ['/Users/sebastienproulx/OneDrive - McGill University/dataBig/C-derived/DecodingHR/fun/y/' subjList{figOption.subj} '/v1wave.mat'];
    warning('need some more coding to show figure without actually processing the thing')
    
    stimPeriod = 12;
    fac = 100;
    Fs = 1;
    n = 120;
    X = ones(n*fac,1);
%     tms = (0:length(X)-1)./Fs;
%     NV = 10;
    minfreq = 1/stimPeriod / 2;
    maxfreq = 1/stimPeriod * 2;
    fb = cwtfilterbank('SignalLength',n*fac,'SamplingFrequency',Fs*fac,...
        'FrequencyLimits',[minfreq maxfreq],...
        'Wavelet','Morse',...
        'TimeBandwidth',4);
    [wave,wave_t] = wavelets(fb);
    [cfs,frq,coi,fb] = cwt(X,'filterbank',fb);
    [~,b] = min(abs(1./frq-stimPeriod));
    1./frq(b)
    thresh = 0.02:0.0001:0.021;
    TimeSupport = nan(size(thresh));
    for i = 1:length(thresh)
        spsi = waveletsupport(fb,thresh(i));
        TimeSupport(i) = spsi(b(1),:).TimeSupport;
    end
    [~,bb] = min(abs(TimeSupport-stimPeriod));
    disp([num2str((1-2*thresh(bb))*100,'%0.2f%%') ' of the energy of the ' num2str(1./frq(b)) 'sec-wavelength wavelet is captured over ' num2str(TimeSupport(bb)) ' seconds around the wavelet''s center']);
    
    figure('WindowStyle','docked');
    hWave = plot(wave_t,imag(wave(b(1),:))); hold on
    plot(xlim,[0 0],':k');
    plot([-0.5 0.5].*stimPeriod,[0 0],'k');
    uistack(hWave,'top');
    spsi = waveletsupport(fb,0.3e-4);
    xlim([spsi(b,:).Begin spsi(b,:).End])
    ax = gca;
    ax.XTick = -24:6:24;
    ax.XGrid = 'on';
    xlabel('time (sec)')
    ylabel('wavelet amplitude (a.u.)')
    return
end

    

for smLevel = {''}
% for smLevel = {'' '_sm1.50' '_sm2.00' '_sm3.00' '_sm4.00' '_sm1.00' '_sm1.25' '_sm1.75' '_sm2.50'}
    for subjInd = 1:length(subjList)
%         try
            
            %% Get data and design
            clearvars -except tstart mask subjInd smLevel subjStimList subjList maskLabel matFun repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion figOption
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
                clearvars -except tstart mask subjInd subjList subjStimList smLevel subj funData_folderIN funData_folderOUT labelDir files label labelList cond sessionLabel data design extraRegr labelList curLabel condCount sessionLabel maskLabel matFun repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion figOption
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
                        nDataPtsToCensor = 0;
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
            
            
            clearvars -except tstart mask subjInd smLevel subjStimList subjList subj funData_folderIN funData_folderOUT labelDir data design extraRegr sessionLabel sessModel stimdur tr maskLabel repo funDir anatDir stimDir inDir outDir noMovement runInd exclude exclusion allTr figOption
            %% GLMdenoise on all sessions (not split)
            % - - - - - - - - - -
            % o x ------------
            splitIn = 1;
            display([subj '; split in ' num2str(splitIn)])
            
            splitPossible = 10;
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
%                 %remove another cycle at the begining, keep the end
%                 splitData{i}(:,:,:,1:12) = [];
%                 splitDesign{i}(1:12) = [];
%                 splitExtraRegr{i}(1:12,:) = [];
                
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
            clearvars -except tstart mask subjInd smLevel subjStimList subjList subj funData_folderIN funData_folderOUT labelDir data design extraRegr sessionLabel sessModel splitDesign splitData splitIn splitExtraRegr splitSessionLabel stimdur tr maskLabel repo funDir anatDir stimDir inDir outDir noMovement runInd splitRunInd exclude exclusion allTr figOption
            
            runInd = 1;
            tmp = permute(splitData{runInd},[4 1 2 3]);
            sz = size(tmp);
            Fs = 1;
            X = mean(tmp(:,:),2)';
            tms = (0:length(X)-1)./Fs;
            NV = 10;
            minfreq = 1/12 / 2;
            maxfreq = 1/12 * 2;
            fb = cwtfilterbank('SignalLength',numel(X),'SamplingFrequency',Fs,...
                'FrequencyLimits',[minfreq maxfreq],...
                'Wavelet','Morse',...
                'TimeBandwidth',4);
            [wave,wave_t] = wavelets(fb);
            [cfs,frq,coi,fb] = cwt(X,'filterbank',fb);
            [~,b] = min(abs(1./frq-12));
            1./frq(b);
            tmsLim = interp1(coi(end/4*3:end),tms(end/4*3:end),frq(b));
            good = tms<tmsLim;
            
            if figOption.subj==subjInd
                figure('WindowStyle','docked');
                ax1 = subplot(3,1,1);
                plot(X); hold on
                ylabel('BOLD signal')
                yyaxis right
                plot(wave_t-wave_t(1)-3,real(wave(b,:)))
                xlim(tms([1 end]));
                ylabel({'Morse wavelet' 'real amplitude'})
                ax2 = subplot(3,1,2);
                hSurf = surface(tms,frq,real(cfs)); hold on
                axis tight
                shading flat
                ylabel('Frequency (Hz)')
                set(gca,'yscale','log')
                hCb = colorbar;
                ylabel(hCb,'real')
                drawnow
                ax2.Position([1 3]) = ax1.Position([1 3]);
                linkaxes([ax1 ax2],'x')
                xlim(tms([1 end]));
                z = ones(size(tms)).*(max(hSurf.ZData(:))*2);
                hP = plot3(tms,ones(size(tms)).*frq(b),z,'k');
                tmsTmp = tms; tmsTmp(coi<min(frq)) = nan;
                coiTmp = coi; coiTmp(coi<min(frq)) = nan;
                plot3(tmsTmp,coiTmp,z,':k')
                ax3 = subplot(3,1,3);
                plot(tms,real(cfs(b,:))); hold on
                ylabel('Signal magnitude')
                yyaxis right
                plot(tms,angle(cfs(b,:))./pi*6); hold on
                ylabel('Signal phase');
                xlabel('Time (s)')
                xlim(tms([1 end]));
                
                %rotate at stim freq
                rho = abs(cfs(b,:));
                theta = angle(cfs(b,:));
                theta = theta - tms.*(2*pi/12);
                [u,v] = pol2cart(theta,rho);
                sig = complex(u,v);
                rho = abs(sig);
                theta = wrapToPi(angle(sig)-angle(mean(sig)));
                [u,v] = pol2cart(theta,rho);
                sig = complex(u,v);
                plot(tms,angle(sig)./pi*6)
                yyaxis left
                plot(tms,real(sig))
                
                hLeg = legend({'real' 'realPhased' 'phase' 'phasePhased'});
                hLeg.Box = 'off';
                if figOption.save
                    saveas(gcf,fullfile(repo,funDir,outDir,[mfilename '.fig']))
                end
            end
            
            
            disp('--- Wavelet Analysis ---')
            ana = 'wave';
            
            splitWave = cell(size(splitData));
            for runInd = 1:length(splitData)
                disp(['Run' num2str(runInd) '/' num2str(length(splitData))])
                tmp = permute(splitData{runInd},[4 1 2 3]); splitData{runInd} = [];
                sig = nan(size(tmp));
                for voxInd = 1:prod(sz(2:end))
                    if all(tmp(:,voxInd)==0)
                    else
                        X = tmp(:,voxInd)';
                        
                        [cfs,~,~,~] = cwt(X,'filterbank',fb);
                        rho = abs(cfs(b,:));
                        theta = angle(cfs(b,:));
                        theta = theta - tms.*(2*pi/12);
                        [u,v] = pol2cart(theta,rho);
                        sig(:,voxInd) = complex(u,v);
                    end
                end
                %store
                splitWave{runInd} = permute(sig,[2 3 4 1]);
            end
            
            splitWave = cell2mat(permute(splitWave,[1 3 4 5 2]));
            results.wave = splitWave; clear splitWave
            results.badEnd = ~good;
            results.badStart = false(size(good)); results.badStart(1:24) = 1;
            
            if noMovement
                if isempty(smLevel{1})
                    outName = fullfile(funData_folderOUT,[maskLabel ana]);
                else
                    outName = fullfile(funData_folderOUT,[maskLabel ana '_' smLevel{1}]);
                end
            else
                if isempty(smLevel{1})
                    outName = fullfile(funData_folderOUT,[maskLabel ana '_move12']);
                else
                    outName = fullfile(funData_folderOUT,[maskLabel ana '_move12' smLevel{1}]);
                end
            end     
            save([outName '.mat'],'results','-v7.3'); clear results % 
            disp([subj '; split in ' num2str(splitIn) ': done'])
            disp(['saved to ' outName])
    end
end
