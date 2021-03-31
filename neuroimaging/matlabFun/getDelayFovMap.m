function [cont,hLine] = getDelayFovMap(d,filterSD,level)
% filterSD in ecc dva
% level between 0 and 1

level = level.*pi;
voxProp = d.voxProp;
voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
[~,U,V,densityXY,X,Y] = pol2surf(voxProp);

ind = true(size(d.sin,1),1);
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

figure('WindowStyle','docked');
[ha, pos] = tight_subplot(2,1,0,0.03,0.03); drawnow
%% Get contour on smoothed map
axes(ha(1))
imagesc(X(1,:),Y(:,1),vecXYsm,[0 pi]); hold on
colormap(flip(cMap,1));
set(gca,'YDir','normal')
hScat = scatter(U,V); hold on
hScat.MarkerEdgeColor = 'none';
hScat.MarkerFaceColor = 'k';
hScat.SizeData = 4^2;
hScat.MarkerFaceAlpha = 0.1;
ax = gca;
ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
xLim = xlim; yLim = ylim;

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
drawnow

%% Plot contour on non-smoothed version
axes(ha(2))
vecXY(isnan(densityXY)) = pi/2;
imagesc(X(1,:),Y(:,1),vecXY,[0 pi]); hold on
colormap(flip(cMap,1));
set(gca,'YDir','normal')
hScat = scatter(U,V); hold on
hScat.MarkerEdgeColor = 'none';
hScat.MarkerFaceColor = 'k';
hScat.SizeData = 4^2;
hScat.MarkerFaceAlpha = 0.1;
ax = gca;
ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
xLim = xlim; yLim = ylim;
for contInd = 1:length(cont)
    plot(cont{contInd}(:,1),cont{contInd}(:,2),'k')
end

hLine = findobj(ha(2).Children,'Type','Line');
[~,b] = sort(cellfun('length',{hLine.XData}),'descend');
hLine = hLine(b);





