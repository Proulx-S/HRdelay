function [indIn,fs,pgon] = processDelayFovContour(d,p,pgonStart,plotFlag)
warning('off','MATLAB:polyshape:repairedBySimplify')
%% Initiate
subjInd = p.subjInd;
sessInd = p.sessInd;

if isfield(p.featSel.fov.empirical,'auto')
    autoSelect1 = 1;
    contInList1 = p.featSel.fov.empirical.auto(1).contInList1;
    autoSelect2 = 1;
    contInList2 = p.featSel.fov.empirical.auto(1).contInList2;
    smList = p.featSel.fov.empirical.auto(1).smList;
    mergeRadiusList = p.featSel.fov.empirical.auto(1).mergeRadiusList;
    marginRadiusList = p.featSel.fov.empirical.auto(1).marginRadiusList;
    
    if length(contInList1)>1
        contInd1 = p.featSel.fov.empirical.contIndList1{sessInd,subjInd};
    else
        contInd1 = contInList1{1};
    end
    if ischar(contInd1)&&strcmp(contInd1,'auto')
        contInd1 = inf;
    end
    
    if length(contInList2)>1
        contInd2 = p.featSel.fov.empirical.contIndList2{sessInd,subjInd};
    else        
        contInd2 = contInList2{1};
    end
    if ischar(contInd2)&&strcmp(contInd2,'auto')
        contInd2 = inf;
    end
    
    if length(smList)>1
        sm = p.featSel.fov.empirical.smList(sessInd,subjInd);
    else
        sm = smList;
    end
    if length(mergeRadiusList)>1
        mergeRadius = p.featSel.fov.empirical.mergeRadiusList(sessInd,subjInd);
    else
        mergeRadius = mergeRadiusList;
    end
    if length(marginRadiusList)>1
        marginRadius = p.featSel.fov.empirical.marginRadiusList(sessInd,subjInd);
    else
        marginRadius = marginRadiusList;
    end
else
    contInd1 = p.featSel.fov.empirical.contIndList1{sessInd,subjInd};
    contInd2 = p.featSel.fov.empirical.contIndList2{sessInd,subjInd};
    sm = p.featSel.fov.empirical.smList(sessInd,subjInd);
    mergeRadius = p.featSel.fov.empirical.mergeRadiusList(sessInd,subjInd);
    marginRadius = p.featSel.fov.empirical.marginRadiusList(sessInd,subjInd);
end

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
    if ~isempty(pgonStart)
        plot(pgonStart,'FaceAlpha',0);
    end
    if ~autoSelect1
        for contInd = 1:pgon.NumRegions
            [x,y] = centroid(pgon,contInd);
            text(x,y,num2str(contInd),'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',12,'Color','m')
        end
    end
    addEccRef(d,p)
    title({['subj' num2str(subjInd) '; sess' num2str(sessInd)] ['sm=' num2str(sm) '; mergeRad=' num2str(mergeRadius) '; marginRad=' num2str(marginRadius)]})
    drawnow
end

if autoSelect1
    if plotFlag>1
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'on');
    else
        f = showDelayFovContour(X,Y,vecXY,cMap,U,V,[],[],'off');
    end
    fs = [fs f];
    [x,y] = addEccRef(d,p);
    Axis = axis;
    pgonRef = polyshape(mat2cell(x,size(x,1),[1 1]),mat2cell(y,size(y,1),[1 1]));
    pgonRef = addboundary(pgonRef,Axis([1 2 2 1 1]),Axis([1 1 2 2 1]));
    pgon = regions(pgon);    
    pgon = pgon(overlaps(pgonRef,pgon));
    [~,b] = sort(area(pgon),'descend');
    pgon = union(pgon(b));
    pgon = simplify(pgon);
    pgon = rmholes(pgon);
    plot(pgon,'FaceAlpha',0)
%     for contInd = 1:pgon.NumRegions
%         [x,y] = centroid(pgon,contInd);
%         text(x,y,num2str(contInd),'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',12,'Color','m')
%     end
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
pgon = union(pgon(b));

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