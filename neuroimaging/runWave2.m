function runWave2(verbose,figOption)
if ~exist('verbose','var')
    verbose = 1;
end

actuallyRun = 1;
if exist('figOption','var') && isfield(figOption,'save')
    saveFig = figOption.save;
else
    saveFig = 0;
end

if ismac
    repo = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repo = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
funDir = 'C-derived\DecodingHR\fun';
inDir = 'y';
outDir = 'y';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repo' 'funDir'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end

% maskLabel = 'v1v2v3';
maskLabel = 'v1';

% subjList = {'03sk'};
% subjStimList = {'sk'};
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
subjStimList = {'jp' 'sk' 'sp' 'bm' 'sb' 'bj'};


if ~actuallyRun
    error('code that')
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


for subjInd = 1:length(subjList)
    disp(['Doing ' subjList{subjInd}])
    filenameIn = fullfile(repo,funDir,inDir,subjList{subjInd},[maskLabel 'SinCos']);
    load([filenameIn '_detrendedData.mat']);
    sz = size(detrendedData{1});
    frqInd = 1;
    cfs = nan([sz(1:3) length(detrendedData) sz(end)]);
    for i = 1:length(detrendedData)
        disp(['Run' num2str(i) '/' num2str(length(detrendedData))])
        [cfsTmp,t,good,wave,wave_t,frq,waveEnergy] = getWaves(detrendedData{i}); detrendedData{i} = [];
        cfs(:,:,:,i,:) = cfsTmp(:,:,:,:,frqInd);
    end
    t = permute(t,[1 2 3 5 4]);
    good = permute(good(:,:,:,:,frqInd),[1 2 3 5 4]);
    wave = permute(wave,[1 2 3 5 4]);
    wave_t = permute(wave_t,[1 2 3 5 4]);
    frq = permute(frq(:,:,:,:,frqInd),[1 2 3 5 4]);
    
    %% Save
    results.wave = cfs; clear cfs
    results.t = t;
    results.good = good;
    results.frq = frq;
    results.waveTD = wave;
    results.waveTD_t = wave_t;
    results.waveTD_eOverPeriod = waveEnergy;
    
    filenameOut = fullfile(repo,funDir,outDir,subjList{subjInd},[maskLabel 'Wave']);
    save([filenameOut '.mat'],'results','-v7.3'); clear results %
    disp(['saved to ' filenameOut '.mat'])
end

function [cfs,t,good,wave,wave_t,frq,waveEnergy] = getWaves(data)

% Prep
data = permute(data,[4 1 2 3]);
sz = size(data);
t = (0:sz(1)-1)';

stimPeriod = 12;
Fs = 1;
n = length(t);
fb = cwtfilterbank('SignalLength',n,'SamplingFrequency',Fs,...
    'FrequencyLimits',[1 2].*1/stimPeriod,...
    'VoicesPerOctave',4,...
    'Wavelet','Morse',...
    'TimeBandwidth',4);
voxInd=1;
[~,frq,coi,~,~] = cwt(data(:,voxInd),'filterbank',fb);
[~,b] = min(abs(frq-1/stimPeriod));
bFund = b(1);
[~,b] = min(abs(frq-1/stimPeriod*2));
bMod = b(1);
frq = frq([bFund bMod],:);
good = true(size(frq));
for i = 1:size(frq,1)
    good(i,t>=interp1(coi(end/4*3:end),t(end/4*3:end),frq(i))) = false;
    good(i,t<=interp1(coi(1:end/4),t(1:end/4),frq(i))) = false;
end

% Wavelet decomposition
cfs = nan([2 sz]);
parfor voxInd = 1:prod(sz(2:end))
    [cfsTmp,~,~,~,~] = cwt(data(:,voxInd),'filterbank',fb);
    cfs(:,:,voxInd) = cfsTmp([bFund bMod],:);
end
cfs = permute(cfs,[3 4 5 2 1]);
frq = permute(frq,[2 3 4 5 1]);
good = permute(good,[3 4 5 2 1]);
t = permute(t,[2 3 4 1 5]);

% Wavelet info
fac = 100;
fb = cwtfilterbank('SignalLength',n*fac,'SamplingFrequency',Fs*fac,...
    'FrequencyLimits',[1 2].*1/stimPeriod,...
    'VoicesPerOctave',4,...
    'Wavelet','Morse',...
    'TimeBandwidth',4);
%wavelet time-domain profile
[wave,wave_t] = wavelets(fb);
wave = permute(wave([bFund bMod],:),[3 4 5 2 1]);
wave_t = permute(wave_t,[3 4 5 2 1]);
spsi = waveletsupport(fb,0.3e-4);
ind = wave_t>=spsi(bFund,:).Begin & wave_t<=spsi(bFund,:).End;
wave_t = wave_t(ind);
wave = wave(ind);
%wavelet time-domain energy at stimulus period
thresh = 0.02:0.0001:0.021;
TimeSupport = nan([1 length(thresh)]);
for i = 1:length(thresh)
    spsi = waveletsupport(fb,thresh(i));
    TimeSupport(:,i) = spsi(bFund,:).TimeSupport;
end
[~,bb] = min(abs(TimeSupport-stimPeriod));
waveEnergy = (1-2*thresh(bb))*100;


