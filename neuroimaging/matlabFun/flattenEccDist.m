function d = flattenEccDist(d,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 1;
end
eccRef = [0.75 7]';
sessInd = 1;

%% Flatten circle dist
nPts = 1000;
ecc = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc{subjInd,sessInd} = d{subjInd,sessInd}.voxProp.ecc;
end
ecc = sort(unique(cat(1,ecc{:})));
cdfecc = sqrt(ecc)/(2*pi);
cdfecc = cdfecc./max(cdfecc);
ecc2 = [0 max(ecc)];
cdfecc2 = [0 1];
eccCirc = interp1(cdfecc2,ecc2,cdfecc);
ecc2eccCirc = smSpline(ecc,eccCirc);

if plotFlag>1
    figure('WindowStyle','docked');
    scatter(ecc,cdfecc); hold on
    scatter(eccCirc,cdfecc); hold on
end
if plotFlag>=1
    subjInd = 1;
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,0)
    title('orig')
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,ecc2eccCirc)
    title('circ')
end

%% Empirically flatten retino dist
ecc = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc{subjInd,sessInd} = ecc2eccCirc(d{subjInd,sessInd}.voxProp.ecc);
end
% Get pdf(ecc) of each subject
if plotFlag>1
    curPlotFlag = 1;
else
    curPlotFlag = 0;
end
eccXY = cell(size(d,1),1);
densityXY = cell(size(d,1),1);
X = cell(size(d,1),1);
Y = cell(size(d,1),1);
okInd = cell(size(d,1),1);
U = cell(size(d,1),1);
V = cell(size(d,1),1);
eccXYspline = cell(size(d,1),1);
densityXYspline = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    [eccXY{subjInd},densityXY{subjInd},X{subjInd},Y{subjInd},okInd{subjInd},U{subjInd},V{subjInd},eccXYspline{subjInd},densityXYspline{subjInd}] = getEccPd(d{subjInd,sessInd},ecc{subjInd,sessInd},curPlotFlag);
    if curPlotFlag
        title(['subj' num2str(subjInd)]);
        figure('WindowStyle','docked');
        hScat = scatter(eccXY{subjInd}(okInd{subjInd}),densityXY{subjInd}(okInd{subjInd}),'ko','filled'); hold on
        alpha(hScat,0.1)
        plot(eccXYspline{subjInd},densityXYspline{subjInd},'r')
        title(['subj' num2str(subjInd)]);
    end
end

% Unifromize pdf(ecc) across subjects
ecc = linspace(min([eccXYspline{:}]),max([eccXYspline{:}]),nPts);
densityecc = nan(size(d,1),nPts);
for subjInd = 1:size(d,1)
    densityecc(subjInd,:) = interp1(eccXYspline{subjInd},densityXYspline{subjInd},ecc);
end
tmpInd = ~all(isnan(densityecc),1);
ecc = ecc(:,tmpInd);
densityecc = densityecc(:,tmpInd);
% Normalize pdf(ecc) according to the overlap across subj
tmpInd = ~any(isnan(densityecc),1);
densityecc = densityecc./sum(densityecc(:,tmpInd),2).*mean(sum(densityecc(:,tmpInd),2));
% Add the average
densityecc = cat(1,densityecc,nanmean(densityecc,1));

if plotFlag>1
    figure('WindowStyle','docked');
    tmpX = cat(3,eccXY{:});
    tmpY = cat(3,densityXY{:});
    tmpInd = cat(3,okInd{:});
    hScat = scatter(tmpX(tmpInd),tmpY(tmpInd),'ko','filled'); hold on
    alpha(hScat,0.03)
    plot(ecc,densityecc(1:end-1,:),':r')
    plot(ecc,densityecc(end,:),'-r')
end

% Get cdf(ecc)
cdfecc = nan(size(d,1),nPts);
for subjInd = 1:size(d,1)+1
    cdfecc(subjInd,~isnan(densityecc(subjInd,:))) = cumsum(densityecc(subjInd,~isnan(densityecc(subjInd,:))),2);
end
cdfecc = cdfecc./max(cdfecc(end,:));
if plotFlag>1
    figure('WindowStyle','docked');
    plot(ecc,cdfecc(1:end-1,:),':r'); hold on
    plot(ecc,cdfecc(end,:),'-r')
end

% Flatten the pdf(ecc) by linearizing the cdf(ecc)
cdfecc = cdfecc(end,:);
ecc2 = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc2{subjInd,sessInd} = d{subjInd,sessInd}.voxProp.ecc;
end
ecc2 = ecc2eccCirc(sort(unique(cat(1,ecc2{:}))));
cdfecc2 = interp1(ecc,cdfecc,ecc2);

ecc = ecc2; clear ecc2
cdfecc = cdfecc2; clear cdfecc2

ecc2 = [0 max(ecc)];
cdfecc2 = [0 1];
eccFlat = interp1(cdfecc2,ecc2,cdfecc);

if plotFlag>1
    figure('WindowStyle','docked');
    scatter(ecc,cdfecc); hold on
    scatter(eccFlat,cdfecc); hold on
end

ecc = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc{subjInd,sessInd} = d{subjInd,sessInd}.voxProp.ecc;
end
ecc = sort(unique(cat(1,ecc{:})));
ecc2eccFlat = smSpline(ecc,eccFlat);

if plotFlag>=1
    subjInd = 1;
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,ecc2eccFlat)
    title('flatten')
end


%% Output
for subjInd = 1:size(d,1)
    d{subjInd,1}.voxProp.ecc2eccFlat = ecc2eccFlat;
    %sess2 is the same as sess1
    d{subjInd,2}.voxProp.ecc2eccFlat = d{subjInd,1}.voxProp.ecc2eccFlat;
end


function [eccXY,densityXY,X,Y,okInd,U,V,eccXYspline,densityXYspline] = getEccPd(d,ecc,plotFlag)
if exist('ecc','var') || ~isempty(ecc)
    d.voxProp.ecc = ecc;
end
if ~exist('plotFlag','var')
    plotFlag = 0;
end
[densityUV,U,V,densityXY,X,Y] = pol2surf(d.voxProp);
if plotFlag
    figure('WindowStyle','docked');
    surf(X,Y,densityXY,'LineStyle','none'); hold on
    scatter3(U,V,densityUV,eps,'w.')
    colormap hot
    alpha(.8)
    set(gca, 'color', 'b');
    ax = gca;
    ax.PlotBoxAspectRatio(1:2) = max(ax.PlotBoxAspectRatio(1:2));
    ax.DataAspectRatio(1:2) = max(ax.DataAspectRatio(1:2));
end
[~,eccXY] = cart2pol(X,Y);
okInd = ~isnan(densityXY);
sm = 0.8;
[fitresult, ~] = fitSpline(eccXY(okInd), densityXY(okInd), sm);
eccXYspline = linspace(min(eccXY(okInd)),max(eccXY(okInd)),100);
densityXYspline = fitresult(eccXYspline);



function [fitresult, gof] = smSpline(curecc, cureccFlat)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( curecc, cureccFlat );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

