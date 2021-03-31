function [hLine,hScat,vecUV] = getDelayFovMap(d,p,ind,filterSD,level,padFac,subjInd)
% filterSD in ecc dva
% level between 0 and 1
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
if ~exist('ind','var') || isempty(ind)
    ind = true(size(d.sin,1),1);
end
if ~exist('padFac','var') || isempty(padFac)
    padFac = 1.2;
end


level = level.*pi;
voxProp = d.voxProp;
if isfield(voxProp,'eccTrans')
    voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
end
voxProp.area = voxProp.area(ind);
voxProp.ecc = voxProp.ecc(ind);
voxProp.pol = voxProp.pol(ind);
voxProp.hemifieldL = voxProp.hemifieldL(ind);
voxProp.hemifieldR = voxProp.hemifieldR(ind);

[~,U,V,densityXY,X,Y] = pol2surf(voxProp,padFac);

vecUV = mean(d.sin(ind,:),2);
[vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
vecUV = abs(angle(vecUV));
F = scatteredInterpolant(U,V,vecUV,'natural','nearest');
vecXY = F(X,Y);
vecXY(isnan(densityXY)) = level;
cMap = redblue(256);
vecSpace = linspace(pi,0,256);
C = interp1(vecSpace,cMap,vecXY,'nearest');




filterSD = filterSD./[mode(diff(X(1,:))) mode(diff(Y(:,1)))];
if length(level)==1
    level = [1 1].*level;
end
if all(filterSD)
    vecXYsm = imgaussfilt(vecXY,filterSD);
else
    vecXYsm = vecXY;
end

if p.figOption.verbose>=3 ...
    || (p.figOption.verbose>=2 && p.figOption.subjInd==subjInd)
    visibilityFlag = 'on';
else
    visibilityFlag = 'off';
end
f = figure('WindowStyle','docked','visible',visibilityFlag);
[ha, ~] = tight_subplot(2,1,0,0.03,0.03); drawnow
%% Get contour on smoothed map
set(f,'CurrentAxes',ha(1))
imagesc(X(1,:),Y(:,1),vecXYsm,[0 pi]); hold on
colormap(flip(cMap,1));
set(ha(1),'YDir','normal')
hScat = scatter(U,V); hold on
hScat.MarkerEdgeColor = 'none';
hScat.MarkerFaceColor = 'k';
hScat.SizeData = 4^2;
hScat.MarkerFaceAlpha = 0.1;
ha(1).DataAspectRatio = ha(1).DataAspectRatio([1 1 3]);
ha(1).PlotBoxAspectRatio = ha(1).PlotBoxAspectRatio([1 1 3]);
xLim = ha(1).XLim; yLim = ha(1).YLim;

[M,c] = contour(X,Y,vecXYsm,level); hold on
c.Visible = 'off';
%read M
cont = cell(0);
contLevel = [];
while ~isempty(M)
    contLevel = [contLevel M(1,1)];
    cont = [cont {M(:,2:1+M(2,1))'}];
    M(:,1:1+M(2,1)) = [];
end
for contInd = 1:length(cont)
    plot(cont{contInd}(:,1),cont{contInd}(:,2),'k')
end
title(['subj' num2str(subjInd)])
drawnow

%% Plot contour on non-smoothed version
vecXY(isnan(densityXY)) = pi/2;
set(f,'CurrentAxes',ha(2))
imagesc(X(1,:),Y(:,1),vecXY,[0 pi]); hold on
colormap(flip(cMap,1));
set(gca,'YDir','normal')
hScat = scatter(U,V); hold on
hScat.MarkerEdgeColor = 'none';
hScat.MarkerFaceColor = 'k';
hScat.SizeData = 4^2;
hScat.MarkerFaceAlpha = 0.1;
ha(2).DataAspectRatio = ha(2).DataAspectRatio([1 1 3]);
ha(2).PlotBoxAspectRatio = ha(2).PlotBoxAspectRatio([1 1 3]);
xLim = ha(2).XLim; yLim = ha(2).YLim;
for contInd = 1:length(cont)
    plot(cont{contInd}(:,1),cont{contInd}(:,2),'k')
end

hLine = findobj(ha(2).Children,'Type','Line');
[~,b] = sort(cellfun('length',{hLine.XData}),'descend');
hLine = hLine(b);




