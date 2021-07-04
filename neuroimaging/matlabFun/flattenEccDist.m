function voxProp = flattenEccDist(d,hemi,p,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 1;
end
eccRef = p.featSel.fov.threshVal';

for sessInd = 1:numel(d)
    d{sessInd} = getOneHemi(d{sessInd},hemi);
end

%% Plot original vox distribution on fov
if plotFlag>=2
    for subjInd = 1:size(d,1)
        plotVoxOnFoV(d{subjInd,p.figOption.sessInd},[],eccRef)
        title(['subj' num2str(subjInd) '; orig'])
    end
elseif plotFlag>=1
    plotVoxOnFoV(d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef)
    title(['subj' num2str(p.figOption.subjInd) '; orig'])
end
%% Flatten vox distribution based on circle area
[ecc2eccFlat_1,eccFlat2ecc_1,~,~,~,~] = getFlatTrans(d,p,'circ',[],plotFlag);
if plotFlag>=2
    for subjInd = 1:size(d,1)
        plotVoxOnFoV(d{subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1{subjInd})
        title(['subj' num2str(subjInd) '; circular flattened'])
    end
elseif plotFlag>=1
    plotVoxOnFoV(d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1{p.figOption.subjInd})
    title(['subj' num2str(p.figOption.subjInd) '; circular flattened'])
end
%% Flatten vox distribution based empirically measured pdf
[ecc2eccFlat_2,eccFlat2ecc_2,~,ecc2eccFlat_1_2,eccFlat2ecc_1_2,~] = getFlatTrans(d,p,'empirical',ecc2eccFlat_1,plotFlag);
if plotFlag>=2
    sessInd = 1;
    for subjInd = 1:size(d,1)
        plotVoxOnFoV(d{subjInd,sessInd},[],eccRef,ecc2eccFlat_1_2{subjInd})
        title(['subj' num2str(subjInd) '; empirically flattened'])
    end
elseif plotFlag>=1
    plotVoxOnFoV(d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1_2{p.figOption.subjInd})
    title(['subj' num2str(p.figOption.subjInd) '; empirically flattened'])
end

% eccXYspline,densityXYspline
% eccGrid,eccGrid_pdf
            
%% Output the ecc transform to d
voxProp = cell(size(d,1),1);
sessInd = 1;
for subjInd = 1:size(d,1)
    voxProp{subjInd} = d{subjInd,sessInd}.voxProp;
    voxProp{subjInd}.eccTrans.toFlat = [ecc2eccFlat_1(subjInd) ecc2eccFlat_2(subjInd) ecc2eccFlat_1_2(subjInd)];
    voxProp{subjInd}.eccTrans.toOrig = [eccFlat2ecc_1(subjInd) eccFlat2ecc_2(subjInd) eccFlat2ecc_1_2(subjInd)];
    voxProp{subjInd}.eccTrans.info = {'circ' 'empirical' 'circ+empirical'};
end

function [ecc2eccFlat,eccFlat2ecc,ecc2eccFlatMean,ecc2eccFlat_fromOrig,eccFlat2ecc_fromOrig,ecc2eccFlatMean_fromOrig] = getFlatTrans(d,p,transFlag,eccTrans,plotFlag)
sessInd = 1;
%% Update ecc with existing flattening transformations
ecc = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    ecc{subjInd} = d{subjInd,sessInd}.voxProp.ecc;
end
if ~isempty(eccTrans)
    eccOrig = cell(size(d,1),1);
    for subjInd = 1:size(d,1)
        eccOrig{subjInd} = ecc{subjInd};
        if isa(eccTrans,'cfit')
            ecc{subjInd} = eccTrans(eccOrig{subjInd});
        else
            ecc{subjInd} = eccTrans{subjInd}(eccOrig{subjInd});
        end
    end
end

switch transFlag
    case 'empirical'
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
        densityXYmed = cell(size(d,1),1);
        eccXYbin = cell(size(d,1),1);
        sessInd = 1;
        sm = 0.99;
        for subjInd = 1:size(d,1)
            [eccXY{subjInd},densityXY{subjInd},X{subjInd},Y{subjInd},okInd{subjInd},U{subjInd},V{subjInd},eccXYspline{subjInd},densityXYspline{subjInd},eccXYbin{subjInd},densityXYmed{subjInd}] = getEccPd(d{subjInd,sessInd},ecc{subjInd,sessInd},sm,plotFlag);
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
        % Replace gaps in pfd of individual subjects with the mean
        for subjInd = 1:size(eccGrid_pdf,1)
            ind = isnan(eccGrid_pdf(subjInd,:));
            eccGrid_pdf(subjInd,ind) = eccGrid_pdfMean(ind);
        end
        % Plot
        if plotFlag>1
            for subjInd = 1:size(d,1)
                if subjInd==1 || plotFlag>2
                    figure('WindowStyle','docked');
                    hScat = scatter(eccXY{subjInd}(okInd{subjInd}),densityXY{subjInd}(okInd{subjInd}),'ko','filled'); hold on
                    alpha(hScat,0.1)
                    plot(eccXYbin{subjInd},densityXYmed{subjInd},'or')
                    plot(eccXYspline{subjInd},densityXYspline{subjInd},'-r')
                    plot(eccGrid,eccGrid_pdf(subjInd,:),'-b')
                    title(['subj' num2str(subjInd)]);
                    legend({'voxels' 'median-filtered' 'spline-smoothed'})
                    xlabel('fov ecc (dva)')
                    ylabel('voxel density on fov')
                end
            end
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
        eccGrid_cdf = nan(size(eccGrid_pdf));
        for subjInd = 1:size(d,1)
            eccGrid_cdf(subjInd,:) = cumsum(eccGrid_pdf(subjInd,:),2);
        end
        eccGrid_cdfMean = cumsum(eccGrid_pdfMean,2);
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
        
%         figure('WindowStyle','docked');
%         subjInd = 6;
%         scatter(ecc{subjInd},ecc_cdf{subjInd}); hold on
%         plot(eccGrid,eccGrid_cdf(subjInd,:))
        
        
        % Interpolate ecc from the empirical cdf(ecc) to a linear cdf(ecc)
        eccFlat = cell(size(d,1),1);
        for subjInd = 1:size(d,1)
            eccLin = [min(ecc{subjInd}) max(ecc{subjInd})];
            ecc_cdfLin = [0 1];
            eccFlat{subjInd} = interp1(ecc_cdfLin,eccLin,ecc_cdf{subjInd});
%             figure('WindowStyle','docked');
%             scatter(ecc{subjInd},ecc_cdf{subjInd}); hold on
%             plot(eccLin,ecc_cdfLin)
%             scatter(eccFlat{subjInd},ecc_cdf{subjInd}); hold on
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


% Define this transformation
[ecc2eccFlat,eccFlat2ecc] = defineTrans(p,ecc,eccFlat,transFlag,plotFlag);
ecc2eccFlatMean = smSpline(cat(1,ecc{:}),eccFlatMean);

if ~isempty(eccTrans)
    % Define the transformation from original ecc to pdf(ecc)-flattened ecc
    [ecc2eccFlat_fromOrig,eccFlat2ecc_fromOrig] = defineTrans(p,eccOrig,eccFlat,transFlag,plotFlag);
    ecc2eccFlatMean_fromOrig = smSpline(cat(1,eccOrig{:}),eccFlatMean);
else
    ecc2eccFlat_fromOrig = ecc2eccFlat;
    eccFlat2ecc_fromOrig = eccFlat2ecc;
    ecc2eccFlatMean_fromOrig = ecc2eccFlatMean;
end

% 
% ecc2eccFlat = cell(size(d,1),1);
% eccFlat2ecc = cell(size(d,1),1);
% eccOrigMax = max(cat(1,eccOrig{:}));
% for subjInd = 1:size(d,1)
%     x = eccOrig{subjInd};
%     y = eccFlat{subjInd};
%     switch transFlag
%         case 'empirical'
%             [a,b] = max(x);
%             y = y./y(b).*a;
%         case 'circ'
%         otherwise
%             error('X')
%     end
%     if plotFlag>1
%         figure('WindowStyle','docked');
%         scatter(x,y); hold on
%     end
%     [x,b] = sort(x);
%     y = y(b);
%     [x,b,~] = unique(x);
%     y = y(b);
%     delta = mean(diff(x(end-10:end)));
%     eccOrigExtra = (x(end)+delta:delta:eccOrigMax)';
%     eccFlatExtra = eccOrigExtra;
%     x = [x; eccOrigExtra];
%     y = [y; eccFlatExtra];
%     ecc2eccFlat{subjInd} = smSpline(x,y);
%     xLim = [min(x) max(x)];
%     x = linspace(xLim(1),xLim(2),1000)';
%     y = ecc2eccFlat{subjInd}(x);
%     eccFlat2ecc{subjInd} = smSpline(y,x);
%     if plotFlag>1
%         plot(x,y); hold on
%         plot(eccFlat2ecc{subjInd}(y),y); hold on
%     end
% end
% 
% 
% % Define the transformation from this transformation to pdf(ecc)-flattened ecc
% ecc2eccFlat_2 = cell(size(d,1),1);
% eccFlat2ecc_2 = cell(size(d,1),1);
% eccMax = max(cat(1,ecc{:}));
% for subjInd = 1:size(d,1)
%     x = ecc{subjInd};
%     y = eccFlat{subjInd};
%     switch transFlag
%         case 'empirical'
%             [a,b] = max(x);
%             y = y./y(b).*a;
%         case 'circ'
%         otherwise
%             error('X')
%     end
%     if plotFlag>1
%         figure('WindowStyle','docked');
%         scatter(x,y); hold on
%     end
%     [x,b] = sort(x);
%     y = y(b);
%     [x,b,~] = unique(x);
%     y = y(b);
%     delta = mean(diff(x(end-10:end)));
%     eccExtra = (x(end)+delta:delta:eccMax)';
%     eccFlatExtra = eccExtra;
%     x = [x; eccExtra];
%     y = [y; eccFlatExtra];
%     ecc2eccFlat_2{subjInd} = smSpline(x,y);
%     xLim = [min(x) max(x)];
%     x = linspace(xLim(1),xLim(2),1000)';
%     y = ecc2eccFlat_2{subjInd}(x);
%     eccFlat2ecc_2{subjInd} = smSpline(y,x);
%     if plotFlag>1
%         plot(x,y); hold on
%         plot(eccFlat2ecc_2{subjInd}(y),y); hold on
%     end
% end
% 
% ecc2eccFlatMean = smSpline(cat(1,eccOrig{:}),eccFlatMean);


function [ecc2eccFlat,eccFlat2ecc] = defineTrans(p,ecc,eccFlat,transFlag,plotFlag)
% Define the transformation from original ecc to pdf(ecc)-flattened ecc
ecc2eccFlat = cell(size(ecc,1),1);
eccFlat2ecc = cell(size(ecc,1),1);
eccMax = max(cat(1,ecc{:}));
for subjInd = 1:size(ecc,1)
    x = ecc{subjInd};
    y = eccFlat{subjInd};
%     switch transFlag
%         case 'empirical'
            [~,b] = min(abs(x-p.featSel.fov.threshVal(2)));
            x1 = x(b);
            if x1<p.featSel.fov.threshVal(2)
                x2 = min(x(x>x1));
            else
                x2 = max(x(x<x1));
            end
            if isempty(x2)
                xRef = p.featSel.fov.threshVal(2);
                param = polyfit(x,y,1);
                yRef = param(1)*p.featSel.fov.threshVal(2)+param(2);
            else
                ind = [find(x1==x,1) find(x2==x,1)];
                xRef = sort(x(ind));
                yRef = sort(y(ind)); clear ind
                yRef = interp1(xRef,yRef,p.featSel.fov.threshVal(2));
                xRef = p.featSel.fov.threshVal(2);
            end
            y = y./yRef.*xRef;            
%         case 'circ'
%         otherwise
%             error('X')
%     end
    if plotFlag>1
        figure('WindowStyle','docked');
        scatter(x,y); hold on
    end
    [x,b] = sort(x);
    y = y(b);
    [x,b,~] = unique(x);
    y = y(b);
    if x(end)~=eccMax
        x = [x; eccMax];
        y = [y; eccMax];
    end
%     y = [y; eccMax.*param(1)+param(2)];
    
%     delta = median(diff(x(round(end*0.90):end)));
%     eccExtra = (x(end)+delta:delta:eccMax)';
%     eccFlatExtra = (y(end)+delta:delta:eccMax)';
% %     eccFlatExtra = eccExtra;
%     x = [x; eccExtra];
%     y = [y; eccFlatExtra];
    ecc2eccFlat{subjInd} = smSpline(x,y);
    xLim = [min(x) max(x)];
    x = linspace(xLim(1),xLim(2),1000)';
    y = ecc2eccFlat{subjInd}(x);
    eccFlat2ecc{subjInd} = smSpline(y,x);
    if plotFlag>1
        plot(x,y); hold on
        plot(eccFlat2ecc{subjInd}(y),y); hold on
    end
end




function [fitresult, gof] = smSpline(curecc, cureccFlat)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( curecc, cureccFlat );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

