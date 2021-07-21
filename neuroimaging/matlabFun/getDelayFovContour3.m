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
% if ~exist('ind','var') || isempty(ind)
%     ind = true(size(d.sin,1),1);
% end

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
