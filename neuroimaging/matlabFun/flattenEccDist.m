function d = flattenEccDist(d,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 1;
end

%% Get pdf(ecc) for each subject
sessInd = 1;
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
pol = cell(size(d,1),1);
U = cell(size(d,1),1);
V = cell(size(d,1),1);
eccXYspline = cell(size(d,1),1);
densityXYspline = cell(size(d,1),1);
for subjInd = 1:size(d,1)
    [eccXY{subjInd},densityXY{subjInd},X{subjInd},Y{subjInd},okInd{subjInd},pol{subjInd},U{subjInd},V{subjInd},eccXYspline{subjInd},densityXYspline{subjInd}] = getEccPd(d{subjInd,sessInd},curPlotFlag);
    if curPlotFlag
        title(['subj' num2str(subjInd)]);
        figure('WindowStyle','docked');
        hScat = scatter(eccXY{subjInd}(okInd{subjInd}),densityXY{subjInd}(okInd{subjInd}),'ko','filled'); hold on
        alpha(hScat,0.1)
        plot(eccXYspline{subjInd},densityXYspline{subjInd},'r')
        title(['subj' num2str(subjInd)]);
    end
end
%% Get the group-averaged pdf(ecc)
nPts = 1000;
ecc = linspace(min([eccXYspline{:}]),max([eccXYspline{:}]),nPts);
densityecc = nan(size(d,1),nPts);
for subjInd = 1:size(d,1)
    densityecc(subjInd,:) = interp1(eccXYspline{subjInd},densityXYspline{subjInd},ecc);
end
if plotFlag>=1
    curPlotFlag = 1;
else
    curPlotFlag = 0;
end
if curPlotFlag
    figure('WindowStyle','docked');
    tmpX = cat(3,eccXY{:});
    tmpY = cat(3,densityXY{:});
    tmpInd = cat(3,okInd{:});
    hScat = scatter(tmpX(tmpInd),tmpY(tmpInd),'ko','filled'); hold on
    alpha(hScat,0.03)
    plot(ecc,densityecc,':r')
    plot(ecc,mean(densityecc,1),'-r')
end
densityecc = median(densityecc,1,'omitnan');
ecc = ecc(~isnan(densityecc));
densityecc = densityecc(~isnan(densityecc));

%% Get the group-averaged cdf(ecc)
cdfecc = cumsum(densityecc);
cdfecc = cdfecc./max(cdfecc);

%% Get the ideal cdf(ecc), one that makes pdf(ecc) flat on a circle
cdfecc2 = linspace(0,1,100);
ecc2 = 2*pi*cdfecc2.^2;
ecc2 = ecc2./max(ecc2).*max(ecc);

%% Put each subject on the group-averaged cdf(ecc) then transform ecc to the ideal cdf(ecc)
if plotFlag>1
    curPlotFlag1 = 1;
else
    curPlotFlag1 = 0;
end
if plotFlag>=1
    curPlotFlag2 = 1;
else
    curPlotFlag2 = 0;
end

for subjInd = 1:size(d,1)
    %group-averaged cdf(ecc)
    curecc = d{subjInd,1}.voxProp.ecc;
    curcdfecc = interp1(ecc,cdfecc,curecc);
    %ideal cdf(ecc)
    cureccFlat = interp1(cdfecc2,ecc2,curcdfecc);
    if curPlotFlag1 && subjInd == 1
        figure('WindowStyle','docked');
        plot(ecc,cdfecc,'k'); hold on
        h1 = plot(curecc,curcdfecc,'ok'); hold on
        plot(ecc2,cdfecc2,'b'); hold on
        h2 = plot(cureccFlat,curcdfecc,'ob'); hold on
        xlabel('ecc (dva)')
        ylabel('density')
        title('flatening circle')
        legend([h1 h2],{'before' 'after'},'location','southeast')
    end
    
    %% Output
    d{subjInd,1}.voxProp.eccFlat.ecc = cureccFlat;
    d{subjInd,1}.voxProp.eccFlat.pol = pol{subjInd};
    d{subjInd,1}.voxProp.eccFlat.U = U{subjInd};
    d{subjInd,1}.voxProp.eccFlat.V = V{subjInd};
    d{subjInd,1}.voxProp.eccXY.X = X{subjInd};
    d{subjInd,1}.voxProp.eccXY.Y = Y{subjInd};
    d{subjInd,1}.voxProp.eccXY.ecc = eccXY{subjInd};
    d{subjInd,1}.voxProp.eccXY.density = densityXY{subjInd};
    d{subjInd,1}.voxProp.eccXY.okInd = okInd{subjInd};
    %sess2 is the same as sess1
    d{subjInd,2}.voxProp.eccFlat = d{subjInd,1}.voxProp.eccFlat;
    d{subjInd,2}.voxProp.eccXY = d{subjInd,1}.voxProp.eccXY;
    
    % Plot
    if curPlotFlag1 || (curPlotFlag2 && subjInd == 1)
        ind = true([size(d{subjInd,sessInd}.sin,1) 1]);
        plotVoxOnFoV(d{subjInd,sessInd},ind,[0.5 8],0)
        title('before flatening')
        plotVoxOnFoV(d{subjInd,sessInd},ind,[],1)
        title('before flatening')
    end
end


function [eccXY,densityXY,X,Y,okInd,pol,U,V,eccXYspline,densityXYspline] = getEccPd(d,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 0;
end
pol = d.voxProp.pol;
pol(d.voxProp.hemifieldL) = d.voxProp.pol(d.voxProp.hemifieldL)./180*pi;
pol(d.voxProp.hemifieldR) = -d.voxProp.pol(d.voxProp.hemifieldR)./180*pi;
pol = wrapToPi(pol+pi/2);
ecc = d.voxProp.ecc;
[U,V] = pol2cart(pol,ecc);
n = 100;
bwFac = 0.009;
delta = (max([U; V]) - min([U; V]))/n;
[X,Y] = meshgrid(linspace(min([U; V])-delta,max([U; V])+delta,n));
[densityXY,~] = ksdensity([U V],[X(:) Y(:)],'Bandwidth',delta*n*bwFac);
densityXY = reshape(densityXY,[n n]);
densityUV = interp2(X,Y,densityXY,U,V);
densityXY(densityXY<(min(densityUV)/10)) = nan;
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

