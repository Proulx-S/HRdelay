function d = getDelayFovContour(d,p,ind)

%% Initialize
U = d.featSel.cont.U;
V = d.featSel.cont.V;
vecUV = d.featSel.cont.vecUV;
X = d.featSel.cont.X;
Y = d.featSel.cont.Y;
outXY = d.featSel.cont.outXY;

sm = p.featSel.fov.empirical.smList(p.sessInd,p.subjInd);
level = p.featSel.fov.empirical.levelList(p.sessInd,p.subjInd);
mergeRadius = p.featSel.fov.empirical.mergeRadiusList(p.sessInd,p.subjInd);
marginRadius = p.featSel.fov.empirical.marginRadiusList(p.sessInd,p.subjInd);
contInd1 = p.featSel.fov.empirical.contIndList1(p.sessInd,p.subjInd);
contInd2 = p.featSel.fov.empirical.contIndList2(p.sessInd,p.subjInd);

level = level.*pi;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
if ~exist('ind','var') || isempty(ind)
    ind = true(size(d.sin,1),1);
end

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
vecXY(outXY|outUV) = level-pi*0.01;

M = contourc(X(1,:),Y(:,1),vecXY,ones(2,1).*level);
cont = cell(0);
while ~isempty(M)
    cont = [cont {M(:,2:1+M(2,1))'}];
    M(:,1:1+M(2,1)) = [];
end

%% Contours to polygons
pgon = polyshape;
for contInd = 1:length(cont)
    pgon = addboundary(pgon,cont{contInd}(:,1),cont{contInd}(:,2));
end
pgon = simplify(pgon);
pgon = rmholes(pgon);
% sort
pgon = regions(pgon);
[~,b] = sort(area(pgon),'descend');
pgon = pgon(b);
pgon2 = polyshape;
for i = 1:length(pgon)
    pgon2 = addboundary(pgon2,pgon(i).Vertices(:,1),pgon(i).Vertices(:,2));
end
pgon = pgon2; clear pgon2


%% Output
d.featSel.cont.vecXY = vecXY;
d.featSel.cont.outUV = outUV;
d.featSel.cont.cMap = cMap;
d.featSel.cont.pgon = pgon;


% contData.cont = cont;
% contData.contLevel = contLevel;
% contData.U = U;
% contData.V = V;
% contData.vecUV = vecUV;
% contData.outUV = outUV;
% contData.X = X;
% contData.Y = Y;
% contData.vecXY = vecXY;
% contData.outXY = outXY;
% contData.cMap = cMap;
