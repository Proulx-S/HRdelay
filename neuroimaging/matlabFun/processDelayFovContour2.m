function [indIn,fs,pgon] = processDelayFovContour2(d,p,sm,mergeRadius,marginRadius,pgonRef,method,plotFlag)
warning('off','MATLAB:polyshape:repairedBySimplify')
%% Initiate
subjInd = p.subjInd;
sessInd = p.sessInd;

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
    plot(pgon,'FaceColor','none');
    [x,y] = addEccRef(d,p);
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
fs = [fs f];
addEccRef(d,p);
% select contours overlapping extended pgonRef
pgonRef = polybuffer(pgonRef,marginRadius);
pgon = regions(pgon);
pgon = union(pgon(overlaps(pgonRef,pgon)));
pgon = simplify(rmholes(pgon));
pgonRef = polybuffer(pgonRef,-marginRadius);
pgon = regions(pgon);
[~,b] = sort(area(pgon),'descend');
pgon = union(pgon(b));
pgon = simplify(pgon);
%         pgon = rmholes(pgon);
plot(pgon,'FaceColor','none')


%% Process contour
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
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
    addEccRef(d,p);
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
    fs = [fs f];
    plot(pgon,'FaceColor','none');
%     plot(pgon,'FaceColor','none');
    addEccRef(d,p);
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

%% Apply contour to feat selection
pgon = regions(pgon);
indIn = false(size(U,1),length(pgon));
for contInd = 1:length(pgon)
    indIn(:,contInd) = inpolygon(U,V,pgon(contInd).Vertices(:,1),pgon(contInd).Vertices(:,2));
end
indIn = any(indIn,2);
pgon = union(pgon);
% plot it
if plotFlag
    if plotFlag>1
        f = showDelayFovContour(X,Y,[],cMap,U(indIn),V(indIn),[],[],'on');
    else
        f = showDelayFovContour(X,Y,[],cMap,U(indIn),V(indIn),[],[],'off');
    end
    fs = [fs f];
    plot(pgon,'FaceColor','none');
    addEccRef(d,p);
    f.Children.XLim = fs(end-1).Children.XLim;
    f.Children.YLim = fs(end-1).Children.YLim;
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end