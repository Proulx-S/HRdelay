function getDelayFovMap(d,filterSD,level,levelColors)
% filterSD in ecc dva
% level between 0 and 1

voxProp = d.voxProp;
voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
[~,U,V,densityXY,X,Y] = pol2surf(voxProp);

ind = true(size(d.sin,1),1);
vecUV = mean(d.sin(ind,:),2);
[vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
vecUV = abs(angle(vecUV));
F = scatteredInterpolant(U,V,vecUV,'natural','nearest');
vecXY = F(X,Y);
%     vecXY(isnan(densityXY)) = nan;
cMap = redblue(256);
vecSpace = linspace(pi,0,256);
C = interp1(vecSpace,cMap,vecXY,'nearest');





figure('WindowStyle','docked');
%     imagesc(X(1,:),Y(:,1),imgaussfilt(vecXY),[0 pi]); hold on
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

filterSD = filterSD./[mode(diff(X(1,:))) mode(diff(Y(:,1)))];
level = level.*pi;
if length(level)==1
    level = [1 1].*level;
end
[M,c] = contour(X,Y,imgaussfilt(vecXY,filterSD),level); hold on
c.Visible = 'off';
%read M
cont = cell(0);
contLevel = [];
while ~isempty(M)
    contLevel = [contLevel M(1,1)];
    cont = [cont {M(:,2:1+M(2,1))'}];
    M(:,1:1+M(2,1)) = [];
end

i = 0;
iCont = 0;
while i<=size(M,2)
    i = i+1;
    %detect new contour
    if ismember(M(1,i),c.LevelList) && round(M(2,i))==M(2,i)
        iCont = iCont + 1;
        cont{iCont} = [];
    else
        cont{iCont} = cat(1,cont{iCont},M(:,i)')
    end
end


M

% c.LineColor = levelColors;
% c.LineColor = 'k';
colormap(gca)

