function d = flattenEccDist(d,p,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 1;
end
eccRef = p.featSel.fov.threshVal';
subjInd = 1;
sessInd = 1;

%% Plot original vox distribution on fov
if plotFlag>=1
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef)
    title('orig')
end
%% Flatten vox distribution based on circle area
[ecc2eccFlat_1,~] = getFlatTrans(d,'circ',[],plotFlag);
if plotFlag>=1
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,ecc2eccFlat_1{subjInd})
    title('circular flattened')
end
%% Flatten vox distribution based empirically measured pdf
[ecc2eccFlat_1_2,~] = getFlatTrans(d,'empiric',ecc2eccFlat_1,plotFlag);
if plotFlag>=1
    plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,ecc2eccFlat_1_2{subjInd})
    title('empirically flattened')
end
%% Output the ecc transform to d
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        d{subjInd,sessInd}.voxProp.eccTrans = ecc2eccFlat_1_2{subjInd};
    end
end

function [ecc2eccFlat,ecc2eccFlatMean] = getFlatTrans(d,transFlag,eccTrans,plotFlag)
sessInd = 1;
%% Update ecc with existing flattening transformations
eccOrig = cell(size(d,1),1);
ecc = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    eccOrig{subjInd} = d{subjInd,sessInd}.voxProp.ecc;
    if ~isempty(eccTrans)
        if isa(eccTrans,'cfit')
            ecc{subjInd} = eccTrans(d{subjInd,sessInd}.voxProp.ecc);
        else
            ecc{subjInd} = eccTrans{subjInd}(d{subjInd,sessInd}.voxProp.ecc);
        end
    else
        ecc{subjInd} = d{subjInd,sessInd}.voxProp.ecc;
    end
end

switch transFlag
    case 'empiric'
        nPts = 1000;
        %% Estimate empirical pdf(ecc)
        eccXY = cell(size(d,1),1);
        densityXY = cell(size(d,1),1);
        X = cell(size(d,1),1);
        Y = cell(size(d,1),1);
        okInd = cell(size(d,1),1);
        U = cell(size(d,1),1);
        V = cell(size(d,1),1);
        eccXYspline = cell(size(d,1),1);
        densityXYspline = cell(size(d,1),1);
        sessInd = 1;
        for subjInd = 1:size(d,1)
            [eccXY{subjInd},densityXY{subjInd},X{subjInd},Y{subjInd},okInd{subjInd},U{subjInd},V{subjInd},eccXYspline{subjInd},densityXYspline{subjInd}] = getEccPd(d{subjInd,sessInd},ecc{subjInd,sessInd},plotFlag);
            if plotFlag>1
                title(['subj' num2str(subjInd)]);
                figure('WindowStyle','docked');
                hScat = scatter(eccXY{subjInd}(okInd{subjInd}),densityXY{subjInd}(okInd{subjInd}),'ko','filled'); hold on
                alpha(hScat,0.1)
                plot(eccXYspline{subjInd},densityXYspline{subjInd},'r')
                title(['subj' num2str(subjInd)]);
            end
        end
        % Interpolate on a grid common across subjects
        eccGrid = linspace(min([eccXYspline{:}]),max([eccXYspline{:}]),nPts);
        eccGrid_pdf = nan(size(d,1),nPts);
        for subjInd = 1:size(d,1)
            eccGrid_pdf(subjInd,:) = interp1(eccXYspline{subjInd},densityXYspline{subjInd},eccGrid);
        end
        tmpInd = ~all(isnan(eccGrid_pdf),1);
        eccGrid = eccGrid(:,tmpInd);
        eccGrid_pdf = eccGrid_pdf(:,tmpInd);
        % Normalize across subjects (according to the region of overlap)
        tmpInd = ~any(isnan(eccGrid_pdf),1);
        eccGrid_pdf = eccGrid_pdf./sum(eccGrid_pdf(:,tmpInd),2).*mean(sum(eccGrid_pdf(:,tmpInd),2));
        % Average across subjects
        eccGrid_pdfMean = nanmean(eccGrid_pdf,1);
        % Plot
        if plotFlag>1
            figure('WindowStyle','docked');
            tmpX = cat(3,eccXY{:});
            tmpY = cat(3,densityXY{:});
            tmpInd = cat(3,okInd{:});
            hScat = scatter(tmpX(tmpInd),tmpY(tmpInd),'ko','filled'); hold on
            alpha(hScat,0.03)
            plot(eccGrid,eccGrid_pdf,':r')
            plot(eccGrid,eccGrid_pdfMean,'-r')
            yLim = ylim; yLim(1) = 0;
            ylim(yLim);
        end
        %% Transform ecc toward a flat pdf(ecc)
        % pdf(ecc) -> cdf(ecc)
        eccGrid_cdf = nan(size(d,1),nPts);
        for subjInd = 1:size(d,1)
            ind = ~isnan(eccGrid_pdf(subjInd,:));
            eccGrid_cdf(subjInd,ind) = cumsum(eccGrid_pdf(subjInd,ind),2);
        end
        ind = ~isnan(eccGrid_pdfMean);
        eccGrid_cdfMean(1,ind) = cumsum(eccGrid_pdfMean(1,ind),2);
        % Normalize cdf(ecc) to max of one based on the group mean (for display
        % purposes)
        eccGrid_cdf = eccGrid_cdf./max(eccGrid_cdfMean);
        eccGrid_cdfMean = eccGrid_cdfMean./max(eccGrid_cdfMean);
        % Plot
        if plotFlag>1
            yyaxis right
            plot(eccGrid,eccGrid_cdf,':b'); hold on
            plot(eccGrid,eccGrid_cdfMean,'-b')
            ax = gca;
            ax.YAxis(1).Color = 'r';
            ax.YAxis(2).Color = 'b';
            yLim = ylim; yLim(1) = 0;
            ylim(yLim);
        end
        % Renormalize cdf(ecc) to max of one on a subj-by-subj basis
        eccGrid_cdf = eccGrid_cdf./max(eccGrid_cdf,[],2);
        eccGrid_cdfMean = eccGrid_cdfMean./max(eccGrid_cdfMean);
        % Interpolate cdf(ecc) at sampled ecc pts
        ecc_cdf = cell(size(d,1),1);
        for subjInd = 1:size(d,1)
            ecc_cdf{subjInd} = interp1(eccGrid,eccGrid_cdf(subjInd,:),ecc{subjInd});
        end
        ecc_cdfMean = interp1(eccGrid,eccGrid_cdfMean,cat(1,ecc{:}));
        % Interpolate ecc from the empirical cdf(ecc) to a linear cdf(ecc)
        eccFlat = cell(size(d,1),1);
        for subjInd = 1:size(d,1)
            eccLin = [min(ecc{subjInd}) max(ecc{subjInd})];
            ecc_cdfLin = [0 1];
            eccFlat{subjInd} = interp1(ecc_cdfLin,eccLin,ecc_cdf{subjInd});
        end
        eccLin = [min(cat(1,ecc{:})) max(cat(1,ecc{:}))];
        ecc_cdfLin = [0 1];
        eccFlatMean = interp1(ecc_cdfLin,eccLin,ecc_cdfMean);
    case 'circ'
        if ~isempty(eccTrans)
            error('does not make sense to perform the circular flattening on flattened ecc')
        end
        eccFlat = cell(size(d,1),1);
        for subjInd = 1:size(d,1)
            ecc_cdf = sqrt(ecc{subjInd})./(2*pi);
            ecc_cdf = ecc_cdf./max(ecc_cdf);
            eccLin = [0 max(ecc{subjInd})];
            ecc_cdfLin = [0 1];
            eccFlat{subjInd} = interp1(ecc_cdfLin,eccLin,ecc_cdf);
        end
        ecc_cdfMean = sqrt(cat(1,ecc{:}))./(2*pi);
        ecc_cdfMean = ecc_cdfMean./max(ecc_cdfMean);
        eccLin = [0 max(cat(1,ecc{:}))];
        ecc_cdfLin = [0 1];
        eccFlatMean = interp1(ecc_cdfLin,eccLin,ecc_cdfMean);
end
if plotFlag>1
    figure('WindowStyle','docked');
    scatter(cat(1,ecc{:}),ecc_cdfMean); hold on
    plot(eccLin,ecc_cdfLin,'k')
    scatter(eccFlatMean,ecc_cdfMean); hold on
end
% Define the transformation from original ecc to pdf(ecc)-flattened ecc
ecc2eccFlat = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc2eccFlat{subjInd} = smSpline(eccOrig{subjInd},eccFlat{subjInd});
end
ecc2eccFlatMean = smSpline(cat(1,eccOrig{:}),eccFlatMean);


function [eccXY,densityXY,X,Y,okInd,U,V,eccXYspline,densityXYspline] = getEccPd(d,ecc,plotFlag)
if exist('ecc','var') || ~isempty(ecc)
    d.voxProp.ecc = ecc;
end
if ~exist('plotFlag','var')
    plotFlag = 0;
end
[densityUV,U,V,densityXY,X,Y] = pol2surf(d.voxProp);
if plotFlag>1
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

