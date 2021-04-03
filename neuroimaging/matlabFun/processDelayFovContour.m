function [indIn,fs] = processDelayFovContour(d,p,plotFlag)
warning('off','MATLAB:polyshape:repairedBySimplify')
%% Initiate
subjInd = p.subjInd;
sessInd = p.sessInd;
contInd1 = p.featSel.fov.empirical.contIndList1{sessInd,subjInd};
contInd2 = p.featSel.fov.empirical.contIndList2{sessInd,subjInd};
sm = p.featSel.fov.empirical.smList(sessInd,subjInd);
mergeRadius = p.featSel.fov.empirical.mergeRadiusList(sessInd,subjInd);
marginRadius = p.featSel.fov.empirical.marginRadiusList(sessInd,subjInd);

X = d.featSel.cont.X;
Y = d.featSel.cont.Y;
vecXY = d.featSel.cont.vecXY;
outXY = d.featSel.cont.outXY;
U = d.featSel.cont.U;
V = d.featSel.cont.V;
vecUV = d.featSel.cont.vecUV;
outUV = d.featSel.cont.outUV;
cMap = d.featSel.cont.cMap;
pgon = d.featSel.cont.pgon;

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
    fs = [fs f];
    plot(pgon,'FaceAlpha',0);
    for contInd = 1:pgon.NumRegions
        [x,y] = centroid(pgon,contInd);
        text(x,y,num2str(contInd),'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',12,'Color','m')
    end
    addEccRef(d,p)
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end


%% Process contour
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
    fs = [fs f];
end
% Remove some contours
if contInd1~=inf
    contInd_rm = 1:pgon.NumRegions;
    contInd_rm = find(~ismember(contInd_rm,contInd1));
    pgon = rmboundary(pgon,contInd_rm);
end
if plotFlag
    plot(pgon,'FaceAlpha',0);
end
% Connect nearby contours by expanding then shrinking them
pgon = polybuffer(pgon,mergeRadius);
pgon = polybuffer(pgon,-mergeRadius);
pgon = simplify(pgon);
pgon = rmholes(pgon);
if plotFlag
    plot(pgon,'FaceAlpha',0);
end
% Expand for conservative boundaries
pgon = polybuffer(pgon,marginRadius);
pgon = simplify(pgon);
pgon = rmholes(pgon);
if plotFlag
    plot(pgon,'FaceAlpha',0);
    addEccRef(d,p)
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

% Sort from larger to smaller area contours
pgon = regions(pgon);
[~,b] = sort(area(pgon),'descend');
pgon = pgon(b);
pgon2 = polyshape;
for i = 1:length(pgon)
    pgon2 = addboundary(pgon2,pgon(i).Vertices(:,1),pgon(i).Vertices(:,2));
end
pgon = pgon2; clear pgon2

% Remove some other contours
if contInd2~=inf
    contInd_rm = 1:pgon.NumRegions;
    contInd_rm = find(~ismember(contInd_rm,contInd1));
    pgon = rmboundary(pgon,contInd_rm);
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
    fs = [fs f];
    plot(pgon,'FaceAlpha',0);
    addEccRef(d,p)
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

%% Apply contour to feat selection
pgon2 = regions(pgon);
cont = {pgon2.Vertices}; clear pgon2
indOut = false(size(d.sin,1),1);
for contInd = 1:length(cont)
    indOut = indOut | inpolygon(U,V,cont{contInd}(:,1),cont{contInd}(:,2));
end
indIn = ~indOut; clear indOut

% plot it
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,[],cMap,U(indIn),V(indIn),[],[],'on');
    else
        f = showDelayFovContour(X,Y,[],cMap,U(indIn),V(indIn),[],[],'off');
    end
    fs = [fs f];
    plot(pgon,'FaceAlpha',0);
    addEccRef(d,p)
    f.Children.XLim = fs(end-1).Children.XLim;
    f.Children.YLim = fs(end-1).Children.YLim;
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end