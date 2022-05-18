function [areaAndFov,cont,voxProp,pEmpirical,fAll,f] = empiricalFov(d,p)

%% Initiate params
pEmpirical.padFac             = 1.2;
pEmpirical.minContPercentArea = 0.05;
pEmpirical.auto(1).smList           = 0.001; % ecc
pEmpirical.auto(1).mergeRadiusList  = 0.70; % ecc
pEmpirical.auto(1).marginRadiusList = 0.40; % ecc
pEmpirical.auto(2).smList           = 0.25; % ecc
pEmpirical.auto(2).mergeRadiusList  = 0.70; % ecc
pEmpirical.auto(2).marginRadiusList = 0.40; % ecc
p.featSel.fov.empirical = pEmpirical;

%% Precompute flattened voxel ecc distribution on fov and delay map
disp('flattening ecc distributions')
voxProp.L = flattenEccDist(d,'L',p);
voxProp.R = flattenEccDist(d,'R',p);
% voxProp.L = flattenEccDist(d,'L',p,1);
% voxProp.R = flattenEccDist(d,'R',p,1);


disp('preparing delay maps')
cont.L = prepareDelayFovContour(d,voxProp.L,p);
cont.R = prepareDelayFovContour(d,voxProp.R,p);

%reorganize voxProp
voxProp2 = cell(size(d));
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        voxProp2{subjInd,sessInd}.L = voxProp.L{subjInd};
        voxProp2{subjInd,sessInd}.R = voxProp.R{subjInd};
    end
    voxProp.L{subjInd} = {};
    voxProp.R{subjInd} = {};
end
voxProp = voxProp2; clear voxProp2


%% Get empirical fov
areaAndFov = cell(size(d));
f.L = cell(size(d));
f.R = cell(size(d));
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        disp([p.meta.subjList{subjInd} ',sess' num2str(sessInd) ': refining fov contour'])
        p.subjInd = subjInd;
        p.sessInd = sessInd;
        hemi = 'L';
        [featValLR.(hemi),featMethodLR.(hemi),featIndInLR.(hemi),featInfoLR.(hemi),f.(hemi){subjInd,sessInd}] = getAreaAndFov(cont.(hemi){subjInd,sessInd},voxProp{subjInd,sessInd}.(hemi),p);
        hemi = 'R';
        [featValLR.(hemi),featMethodLR.(hemi),featIndInLR.(hemi),featInfoLR.(hemi),f.(hemi){subjInd,sessInd}] = getAreaAndFov(cont.(hemi){subjInd,sessInd},voxProp{subjInd,sessInd}.(hemi),p);
        
        % Combine hemifields
        featVal = nan(size(voxProp{subjInd,sessInd}.(hemi).hemifield));
        featIndIn = nan(size(voxProp{subjInd,sessInd}.(hemi).hemifield));
        hemi = 'L';
        featVal(voxProp{subjInd,sessInd}.(hemi).hemifield) = featValLR.(hemi);
        featIndIn(voxProp{subjInd,sessInd}.(hemi).hemifield) = featIndInLR.(hemi);
        hemi = 'R';
        featVal(voxProp{subjInd,sessInd}.(hemi).hemifield) = featValLR.(hemi);
        featIndIn(voxProp{subjInd,sessInd}.(hemi).hemifield) = featIndInLR.(hemi);
        featMethod = featMethodLR.L;
        featInfo = featInfoLR.L;
        
        areaAndFov{subjInd,sessInd}.featVal = featVal; clear featVal featValLR;
        areaAndFov{subjInd,sessInd}.featIndIn = featIndIn; clear featIndIn featIndInLR;
        areaAndFov{subjInd,sessInd}.featMethod = featMethod; clear featMethod featMethodLR;
        areaAndFov{subjInd,sessInd}.featInfo = featInfo; clear featInfo featInfoLR;
    end
end

fIndList = [1 2 4 7 9 10];
supTitleList = {'removing small islands' '1st contours' '1st contours processing' '2nd contours' '2nd contours processing' 'Final contours'};
fAll = cell(size(fIndList));
for i = 1:length(fIndList)
    fInd = fIndList(i);
%     if fIndList(i)==10 || p.figOption.verbose>=1
%         fAll{i} = figure('WindowStyle','docked','visible','on');
%     else
        fAll{i} = figure('WindowStyle','docked','visible','off');
%     end
    [ha, pos] = tight_subplot(size(d,2), size(d,1)*2, 0, 0.1, 0); delete(ha);
    for subjInd = 1:size(d,1)
        for sessInd = 1:size(d,2)
            hemi = 'L';
            ax.(hemi) = copyobj(f.(hemi){subjInd,sessInd}(fInd).Children,fAll{i});
            ax.(hemi).DataAspectRatioMode = 'auto';
            ax.(hemi).PlotBoxAspectRatioMode = 'auto';
            ax.(hemi).Position = pos{(sessInd-1)*(size(d,1)*2)+(subjInd*2-1)};
            ax.(hemi).Colormap = f.(hemi){subjInd,sessInd}(fInd).Children.Colormap;
            drawnow
            
            hemi = 'R';
            ax.(hemi) = copyobj(f.(hemi){subjInd,sessInd}(fInd).Children,fAll{i});
            ax.(hemi).DataAspectRatioMode = 'auto';
            ax.(hemi).PlotBoxAspectRatioMode = 'auto';
            ax.(hemi).Position = pos{(sessInd-1)*(size(d,1)*2)+(subjInd*2-1)+1};
            ax.(hemi).Colormap = f.(hemi){subjInd,sessInd}(fInd).Children.Colormap;
            drawnow
            
            yLim = [-1 1].*max(abs([ax.L.YLim ax.R.YLim ax.L.XLim(1) ax.R.XLim(2)]));
            xLim = yLim(2);
            ax.L.YLim = yLim;
            ax.R.YLim = yLim;
            ax.L.XLim = [-xLim 0];
            ax.R.XLim = [0 xLim];
            drawnow
            
            ax.L.PlotBoxAspectRatio = [0.5 1 1];
            ax.R.PlotBoxAspectRatio = [0.5 1 1];
            
            ax.L.YAxis.Visible = 'off';
            ax.R.YAxis.Visible = 'off';
        end
    end
    sgtitle(supTitleList{i})
end



function voxProp = flattenEccDist(d,hemi,p,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = p.figOption.verbose;
end
eccRef = p.fov.eccLim';

for sessInd = 1:numel(d)
    d{sessInd} = getOneHemi(d{sessInd},hemi);
end

%% Plot original vox distribution on fov
% if plotFlag>=2
%     for subjInd = 1:size(d,1)
%         plotVoxOnFoV(p,d{subjInd,p.figOption.sessInd},[],eccRef)
%         title(['subj' num2str(subjInd) '; orig'])
%     end
if plotFlag>=1
    plotVoxOnFoV(p,d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef)
    title(['subj' num2str(p.figOption.subjInd) '; orig'])
end
%% Flatten vox distribution based on circle area
[ecc2eccFlat_1,eccFlat2ecc_1,~,~,~,~] = getFlatTrans(d,p,'circ',[],plotFlag);
% if plotFlag>=2
%     for subjInd = 1:size(d,1)
%         plotVoxOnFoV(p,d{subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1{subjInd})
%         title(['subj' num2str(subjInd) '; circular flattened'])
%     end
if plotFlag>=1
    plotVoxOnFoV(p,d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1{p.figOption.subjInd})
    title(['subj' num2str(p.figOption.subjInd) '; circular flattened'])
end
%% Flatten vox distribution based empirically measured pdf
[ecc2eccFlat_2,eccFlat2ecc_2,~,ecc2eccFlat_1_2,eccFlat2ecc_1_2,~] = getFlatTrans(d,p,'empirical',ecc2eccFlat_1,plotFlag);
% if plotFlag>=2
%     sessInd = 1;
%     for subjInd = 1:size(d,1)
%         plotVoxOnFoV(p,d{subjInd,sessInd},[],eccRef,ecc2eccFlat_1_2{subjInd})
%         title(['subj' num2str(subjInd) '; empirically flattened'])
%     end
if plotFlag>=1
    plotVoxOnFoV(p,d{p.figOption.subjInd,p.figOption.sessInd},[],eccRef,ecc2eccFlat_1_2{p.figOption.subjInd})
    title(['subj' num2str(p.figOption.subjInd) '; empirically flattened'])
end

            
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


function [ecc2eccFlat,eccFlat2ecc] = defineTrans(p,ecc,eccFlat,transFlag,plotFlag)
% Define the transformation from original ecc to pdf(ecc)-flattened ecc
ecc2eccFlat = cell(size(ecc,1),1);
eccFlat2ecc = cell(size(ecc,1),1);
eccMax = max(cat(1,ecc{:}));
for subjInd = 1:size(ecc,1)
    x = ecc{subjInd};
    y = eccFlat{subjInd};
    [~,b] = min(abs(x-p.fov.eccLim(2)));
    x1 = x(b);
    if x1<p.fov.eccLim(2)
        x2 = min(x(x>x1));
    else
        x2 = max(x(x<x1));
    end
    if isempty(x2)
        xRef = p.fov.eccLim(2);
        param = polyfit(x,y,1);
        yRef = param(1)*p.fov.eccLim(2)+param(2);
    else
        ind = [find(x1==x,1) find(x2==x,1)];
        xRef = sort(x(ind));
        yRef = sort(y(ind)); clear ind
        yRef = interp1(xRef,yRef,p.fov.eccLim(2));
        xRef = p.fov.eccLim(2);
    end
    y = y./yRef.*xRef;
    
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


function d = getOneHemi(d,hemi)

ind = d.voxProp.(['hemifield' hemi]);

fieldList = fields(d);
for fieldInd = 1:length(fieldList)
    sz = size(d.(fieldList{fieldInd}));
    if any(sz==length(ind))
        sz(sz==length(ind)) = nnz(ind);
        tmp = nan(sz);
        sz(sz==nnz(ind)) = 1;
        tmp(:) = d.(fieldList{fieldInd})(repmat(ind,sz));
        d.(fieldList{fieldInd}) = tmp; clear tmp sz
    end
end

fieldList = fields(d.voxProp);
for fieldInd = 1:length(fieldList)
    sz = size(d.voxProp.(fieldList{fieldInd}));
    if any(sz==length(ind))
        sz(sz==length(ind)) = nnz(ind);
        tmp = nan(sz);
        sz(sz==nnz(ind)) = 1;
        tmp(:) = d.voxProp.(fieldList{fieldInd})(repmat(ind,sz));
        d.voxProp.(fieldList{fieldInd}) = tmp; clear tmp sz
    end
end

d.voxProp.hemifieldL = logical(d.voxProp.hemifieldL);
d.voxProp.hemifieldR = logical(d.voxProp.hemifieldR);
d.voxProp.hemifield = ind;


function plotVoxOnFoV(p,d,ind,eccRef,flatFlag)
if ~exist('ind','var') || isempty(ind)
    ind = true(size(d.sin,1),1);
end
if ~exist('eccRef','var')
    eccRef = [];
end
if ~exist('flatFlag','var')
    flatFlag = 0;
elseif isa(flatFlag,'cfit')
    flatTrans = flatFlag;
    flatFlag = 1;
elseif flatFlag
    flatTrans = d.voxProp.eccTrans;
end

vec = mean(d.sin(ind,:),2);
[vec,~] = polarSpaceNormalization(vec,'cartRoi');
%         plotDensity(vec)
vec = abs(angle(vec));

cMap = redblue(256);
vecSpace = linspace(pi,0,256);
c = interp1(vecSpace,cMap,vec,'nearest');

figure('WindowStyle','docked');
ecc = d.voxProp.ecc(ind);
if flatFlag
    ecc = flatTrans(ecc);
end
pol = d.voxProp.pol(ind)./180*pi;
hemiL = logical(d.voxProp.hemifieldL(ind));
if any(hemiL)
    polarscatter(pol(hemiL),ecc(hemiL),eps,c(hemiL,:),'.'); hold on
end
hemiR = logical(d.voxProp.hemifieldR(ind));
if any(hemiR)
    polarscatter(-pol(hemiR),ecc(hemiR),eps,c(hemiR,:),'.'); hold on
end
ax = gca;
ax.ThetaZeroLocation = 'top';


if ~isempty(eccRef)
    th = repmat(linspace(0,2*pi,50),[2 1])';
    if flatFlag
        r = flatTrans(eccRef)';
    else
        r = eccRef';
    end
    polarplot(th,r+zeros(size(th)),'k');
end

if flatFlag
    ax = gca;
    ticks = 0:1:eccRef(2);
    ax.RTick = flatTrans(ticks);
    tickLabels = cellstr(num2str(ticks'));
    tickLabels(9:end) = {''};
    ax.RTickLabel = tickLabels;
end


function [eccXY,densityXY,X,Y,okInd,U,V,eccXYspline,densityXYspline,eccXYbin,densityXYmed] = getEccPd(d,ecc,sm,plotFlag)
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
if exist('sm','var') && ~isempty(sm)
    sm2 = sm;
end
sm = 0.8;
%median filter
n = 100;
if ~any(d.voxProp.hemifieldL) || ~any(d.voxProp.hemifieldR)
    n = n/2;
end
[eccXYbin,densityXYmed] = medianFilter(eccXY(okInd),densityXY(okInd),n);

%smothing spline
% [fitresult, ~] = fitSpline(eccXY(okInd), densityXY(okInd), sm);
[fitresult, ~] = fitSpline(eccXYbin, densityXYmed, sm2);
eccXYspline = linspace(min(eccXY(okInd)),max(eccXY(okInd)),100);
densityXYspline = fitresult(eccXYspline);
densityXYspline(densityXYspline<0) = nan;


function [fitresult, gof] = fitSpline(rho, densityXY, sm, plotFlag)
warning('off','curvefit:prepareFittingData:removingNaNAndInf');
if ~exist('plotFlag','var')
    plotFlag = 0;
end

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( rho, densityXY );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
if exist('sm','var') && ~isempty(sm)
    opts.SmoothingParam = sm;
else
    opts.SmoothingParam = 0.8;
end

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if plotFlag
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'densityXY vs. rho', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'rho', 'Interpreter', 'none' );
    ylabel( 'densityXY', 'Interpreter', 'none' );
    grid on
end


function [densityUV,U,V,densityXY,X,Y] = pol2surf(voxProp,padFac,n)
if ~exist('padFac','var') || isempty(padFac)
    padFac = 1.2;
end
if ~exist('n','var') || isempty(n)
    n = 2^8;
end
ecc = voxProp.ecc;
pol = voxProp.pol;
if any(voxProp.hemifieldL)
    pol(voxProp.hemifieldL) = voxProp.pol(voxProp.hemifieldL)./180*pi;
end
if any(voxProp.hemifieldR)
    pol(voxProp.hemifieldR) = -voxProp.pol(voxProp.hemifieldR)./180*pi;
end
pol = wrapToPi(pol+pi/2);
[U,V] = pol2cart(pol,ecc);

bwFac = 0.009;
delta = range([U; V])*padFac/n;
[X,Y] = meshgrid(linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),n);
[densityXY,~] = ksdensity([U V],[X(:) Y(:)],'Bandwidth',delta*n*bwFac);
densityXY = reshape(densityXY,[n n]);
densityUV = interp2(X,Y,densityXY,U,V);
densityXY(densityXY<(min(densityUV)/10)) = nan;


function [xbin,ymed] = medianFilter(x,y,n)

[~,edges,bin] = histcounts(x,n);
xbin = edges(1:end-1)+diff(edges([1 end]))/n/2;
ymed = nan(size(xbin));
for binInd = 1:n
    ymed(binInd) = median(y(bin==binInd));
end



function cont = prepareDelayFovContour(d,voxProp,p)
padFac = p.featSel.fov.empirical.padFac;
cont = cell(size(d));
for subjInd = 1:size(d,1)
    curVoxProp = voxProp{subjInd};
    % Apply precomputed flattening
    if isfield(curVoxProp,'eccTrans')
        curVoxProp.ecc = curVoxProp.eccTrans.toFlat{end}(curVoxProp.ecc);
    end
    for sessInd = 1:size(d,2)
        vecUV = mean(d{subjInd,sessInd}.sin(curVoxProp.hemifield,:),2);
        
        % Prepare surface
        [~,U,V,densityXY,X,Y] = pol2surf(curVoxProp,padFac);
        outXY = isnan(densityXY); clear densityXY
        [vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
        vecUV = abs(angle(vecUV));
        
        % Output
        cont{subjInd,sessInd}.U = U;
        cont{subjInd,sessInd}.V = V;
        cont{subjInd,sessInd}.vecUV = vecUV;
        cont{subjInd,sessInd}.X = X;
        cont{subjInd,sessInd}.Y = Y;
        cont{subjInd,sessInd}.outXY = outXY;
    end
end


function [featVal,featMethod,featIndIn,featInfo,f] = getAreaAndFov(cont,voxProp,p)

warning('off','MATLAB:polyshape:repairedBySimplify')
if p.figOption.verbose==0
    visibleFlag = 0;
elseif p.figOption.verbose>=1
    visibleFlag = 1;
end

%% Precompute stats on response vector (random effect)
ind = true(size(cont.U,1),1);

%% Voxels representing stimulus fov
curInfo1 = {'fov'};
minContPercentArea = p.featSel.fov.empirical.minContPercentArea;
M = contourc(cont.X(1,:),cont.Y(:,1),double(cont.outXY),ones(2,1).*0.5);
pgon = polyshape;
while ~isempty(M)
    pgon = addboundary(pgon,M(1,2:1+M(2,1)),M(2,2:1+M(2,1)));
    M(:,1:1+M(2,1)) = [];
end
pgonOrig = pgon;

pgon = regions(pgon);
areas = area(pgon);
pgon = pgon(areas/sum(areas) > minContPercentArea);

tmpInd = false([size(cont.X) length(pgon)]);
tmpInd2 = false([size(cont.U,1) length(pgon)]);
for i = 1:length(pgon)
    Vertices = pgon(i).Vertices;
    tmpInd(:,:,i) = inpolygon(cont.X,cont.Y,Vertices(:,1),Vertices(:,2));
    tmpInd2(:,i) = inpolygon(cont.U,cont.V,Vertices(:,1),Vertices(:,2));
end
tmpInd = any(tmpInd,3);
tmpInd2 = any(tmpInd2,2);
if visibleFlag
    visibleFlag2 = 2;
else
    visibleFlag2 = 1;
end
if visibleFlag
    f0 = figure('WindowStyle','docked','visible','on');
else
    f0 = figure('WindowStyle','docked','visible','off');
end
imagesc(cont.X(1,:),cont.Y(:,1),~tmpInd); hold on
set(gca,'YDir','normal'); colormap autumn
plot(pgonOrig,'FaceColor','none')
set(gca,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1]);

cont.outXY = ~tmpInd;
prevInd = ind&tmpInd2;

sm           = p.featSel.fov.empirical.auto(1).smList;
mergeRadius  = p.featSel.fov.empirical.auto(1).mergeRadiusList;
marginRadius = p.featSel.fov.empirical.auto(1).marginRadiusList;
cont = getDelayFovContour3(cont,sm,prevInd);
featVal = cont.vecUV;
[~,f1,pgon] = processDelayFovContour3(cont,voxProp,p,prevInd,sm,mergeRadius,marginRadius,[],'do not add pgonRef',visibleFlag2);

if p.featSel.fov.empirical.auto(2).smList~=sm
    sm = p.featSel.fov.empirical.auto(2).smList;
    cont = getDelayFovContour3(cont,sm,prevInd);
end
mergeRadius  = p.featSel.fov.empirical.auto(2).mergeRadiusList;
marginRadius = p.featSel.fov.empirical.auto(2).marginRadiusList;
if ~isempty(pgon)
    [curIndIn,f2,~] = processDelayFovContour3(cont,voxProp,p,prevInd,sm,mergeRadius,marginRadius,pgon,'add pgonRef',visibleFlag2);
else
    curIndIn = true(size(prevInd));
    f2 = gobjects([1 5]);
end
f = [f0 f1 f2];
featVal;
featMethod = strjoin([{'fov'} {'empirical'}],': ');
featIndIn = curIndIn;
featInfo = 'vox X featSel X condPair';



function cont = getDelayFovContour3(cont,sm,ind)

%% Initialize
U = cont.U;
V = cont.V;
vecUV = cont.vecUV;
X = cont.X;
Y = cont.Y;
outXY = cont.outXY;

level = 0.5;
level = level.*pi;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

%% Interpolate delay to surface
F = scatteredInterpolant(U(ind),V(ind),ones(size(U(ind))).*pi/2,'nearest','none');
outUV = isnan(F(X,Y));

F = scatteredInterpolant(U(ind),V(ind),vecUV(ind),'natural','none');
vecXY = F(X,Y);
cMap = redblue(256);
% figure('WindowStyle','docked');
% hIm = imagesc(X(1,:),Y(:,1),vecXY); hold on
% scatter(U,V,'.k')

vecXY(outXY|outUV) = nan;
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),vecXY,hIm.Parent.CLim); hold on
% scatter(U,V,'.k')

%% Smooth
sigma = sm./mode(diff(X(1,:)));
kernel = fspecial('gauss',2*ceil(2*sigma.*[1 1])+1,sigma);
vecXY = nanconv(vecXY,kernel,'nanout');

%% Get contours
vecXY(outXY|outUV) = level-pi*0.01; %just to out-of-FOV regions behave

M = contourc(X(1,:),Y(:,1),vecXY,ones(2,1).*level);
if isempty(M)
    pgon = cont.pgon;
else
    cont2 = cell(0);
    while ~isempty(M)
        cont2 = [cont2 {M(:,2:1+M(2,1))'}];
        M(:,1:1+M(2,1)) = [];
    end
    
    %% Contours to polygons
    pgon = polyshape;
    for contInd = 1:length(cont2)
        pgon = addboundary(pgon,cont2{contInd}(:,1),cont2{contInd}(:,2));
    end
    pgon = simplify(pgon);
    pgon = rmholes(pgon);
    % sort
    pgon = regions(pgon);
    [~,b] = sort(area(pgon),'descend');
    pgon = union(pgon(b));
end

%% Output
cont.vecXY = vecXY;
cont.outUV = outUV;
cont.cMap = cMap;
cont.pgon = pgon;
cont.ind = ind;
cont.F = F;



function [indIn,fs,pgon] = processDelayFovContour3(cont,voxProp,p,prevInd,sm,mergeRadius,marginRadius,pgonRef,method,plotFlag)
warning('off','MATLAB:polyshape:repairedBySimplify')
%% Initiate
subjInd = p.subjInd;
sessInd = p.sessInd;

X = cont.X;
Y = cont.Y;
vecXY = cont.vecXY;
outXY = cont.outXY;
U = cont.U(prevInd);
V = cont.V(prevInd);
vecUV = cont.vecUV(prevInd);
outUV = cont.outUV;
cMap = cont.cMap;
pgon = cont.pgon;

fs = [];


%% Show original contours
% F = scatteredInterpolant(contData.U,contData.V,contData.vecUV,interp,extrap);
% vecXY = F(contData.X,contData.Y);
vecXY(outXY|outUV) = pi/2;
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
    f.UserData.F = cont.F;
%     f.UserData = cont;
%     f.UserData.F = cont.F;
%     f.UserData.U = cont.U;
%     f.UserData.V = cont.V;
%     f.UserData.vecUV = cont.vecUV;
%     f.UserData.outUV = cont.outUV;
%     f.UserData.X = cont.X;
%     f.UserData.Y = cont.Y;
%     f.UserData.vecXY = cont.vecXY;
%     f.UserData.outXY = cont.vecXY;
    f.UserData.ind = cont.ind;
    fs = [fs f];
    plot(pgon,'FaceColor','none');
    [x,y] = addEccRef(voxProp,p);
    if isempty(pgonRef)
        Axis = axis;
        pgonRef = polyshape(mat2cell(x,size(x,1),[1 1]),mat2cell(y,size(y,1),[1 1]));
        pgonRef = addboundary(pgonRef,Axis([1 2 2 1 1]),Axis([1 1 2 2 1]));
    end
    plot(pgonRef,'FaceColor','y');
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

%% Show automatically slected contours
if plotFlag>1
    f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
else
    f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
end
f.UserData.F = cont.F;
% f.UserData.U = cont.U;
% f.UserData.V = cont.V;
% f.UserData.vecUV = cont.vecUV;
% f.UserData.X = cont.X;
% f.UserData.Y = cont.Y;
% f.UserData.vecXY = cont.vecXY;
f.UserData.ind = cont.ind;
fs = [fs f];
addEccRef(voxProp,p);
% select contours overlapping extended pgonRef
pgonRef = polybuffer(pgonRef,marginRadius);
try
    pgon = regions(pgon);
catch
    keyboard
end
pgon = union(pgon(overlaps(pgonRef,pgon)));
pgon = simplify(rmholes(pgon));
pgonRef = polybuffer(pgonRef,-marginRadius);
if isempty(pgon)
    % Stop here if no region is selected for exclusion
    indIn = true(length(prevInd),1);
    fs = [fs f];
    fs = [fs f];
    fs = [fs f];
    return
else
    pgon = regions(pgon);
    [~,b] = sort(area(pgon),'descend');
    pgon = union(pgon(b));
    pgon = simplify(pgon);
    %         pgon = rmholes(pgon);
    plot(pgon,'FaceColor','none')
end


%% Process contour
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
    f.UserData.F = cont.F;
%     f.UserData.U = cont.U;
%     f.UserData.V = cont.V;
%     f.UserData.vecUV = cont.vecUV;
%     f.UserData.X = cont.X;
%     f.UserData.Y = cont.Y;
%     f.UserData.vecXY = cont.vecXY;
    f.UserData.ind = cont.ind;
    fs = [fs f];
    plot(pgon,'FaceColor','none');
end
% Connect nearby contours by expanding then shrinking them
pgon = polybuffer(pgon,mergeRadius);
pgon = polybuffer(pgon,-mergeRadius);
pgon = simplify(pgon);
pgon = regions(pgon);
[~,b] = sort(area(pgon),'descend');
pgon = union(pgon(b));
if plotFlag
    plot(pgon,'FaceColor','none');
end
% % Select the biggest one of the merged center coutours
% pgonRef2 = regions(pgonRef);
% [~,b] = min(area(pgonRef2));
% pgonRef2 = polybuffer(pgonRef2(b),marginRadius);
% pgon = regions(pgon);
% overInd = find(overlaps(pgonRef2,pgon));
% pgon(overInd(2:end)) = [];
% pgon = union(pgon);

% Expand for conservative boundaries
pgon = polybuffer(pgon,marginRadius);
pgon = simplify(pgon);
switch method
    case 'add pgonRef'
        % add pgonRef to contours
        pgon = simplify(union(pgon,pgonRef));
        pgon = regions(pgon);
        [~,b] = sort(area(pgon),'descend');
        pgon = union(pgon(b));
    case 'do not add pgonRef'
    otherwise
end
if plotFlag
    plot(pgon,'FaceColor','none');
    addEccRef(voxProp,p);
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

%% Show processed contours
F = scatteredInterpolant(U,V,vecUV,'nearest','none');
vecXY = F(X,Y);
vecXY(outXY|outUV) = pi/2;
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
    f.UserData.F = F;
%     f.UserData.U = cont.U;
%     f.UserData.V = cont.V;
%     f.UserData.vecUV = cont.vecUV;
%     f.UserData.X = cont.X;
%     f.UserData.Y = cont.Y;
%     f.UserData.vecXY = cont.vecXY;
    f.UserData.ind = prevInd;
    fs = [fs f];
    plot(pgon,'FaceColor','none');
%     plot(pgon,'FaceColor','none');
    addEccRef(voxProp,p);
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

%% Apply contour to feat selection
pgon = regions(pgon);
indIn = false(length(prevInd),length(pgon));
% indIn = false(size(U,1),length(pgon));
for contInd = 1:length(pgon)
    indIn(prevInd,contInd) = inpolygon(U,V,pgon(contInd).Vertices(:,1),pgon(contInd).Vertices(:,2));
end
indIn = ~any(indIn,2);
pgon = union(pgon);
% plot it
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,[],cMap,U(indIn(prevInd)),V(indIn(prevInd)),[],[],'on');
    else
        f = showDelayFovContour(X,Y,[],cMap,U(indIn(prevInd)),V(indIn(prevInd)),[],[],'off');
    end
    fs = [fs f];
    plot(pgon,'FaceColor','none');
    addEccRef(voxProp,p);
    f.Children.XLim = fs(end-1).Children.XLim;
    f.Children.YLim = fs(end-1).Children.YLim;
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end


function [x,y] = addEccRef(voxProp,p)
theta = linspace(-pi,pi,100);
eccRef = voxProp.eccTrans.toFlat{end}(p.fov.eccLim);
x = nan(size(theta,2),length(eccRef));
y = nan(size(theta,2),length(eccRef));
for i = 1:length(eccRef)
    [x(:,i),y(:,i)] = pol2cart(theta,eccRef(i));
%     plot(x,y,'--','Color',[1 1 1].*0.4)
    plot(x(:,i),y(:,i),'--k')
end