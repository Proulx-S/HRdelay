function contData = getDelayFovContour(d,p,ind,filterSD,level,padFac,subjInd,interp,extrap)
% filterSD in ecc dva
% level between 0 and 1
level = level.*pi;

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
if ~exist('ind','var') || isempty(ind)
    ind = true(size(d.sin,1),1);
end
if ~exist('padFac','var') || isempty(padFac)
    padFac = 1.2;
end

voxProp = d.voxProp;
if isfield(voxProp,'eccTrans')
    voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
end
% voxProp.area = voxProp.area(ind);
% voxProp.ecc = voxProp.ecc(ind);
% voxProp.pol = voxProp.pol(ind);
% voxProp.hemifieldL = voxProp.hemifieldL(ind);
% voxProp.hemifieldR = voxProp.hemifieldR(ind);

[~,U,V,densityXY,X,Y] = pol2surf(voxProp,padFac);
outXY = isnan(densityXY); clear densityXY
F = scatteredInterpolant(U(ind),V(ind),ones(size(U(ind))).*pi/2,'nearest','none');
outUV = isnan(F(X,Y));
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),outUV); hold on
% scatter(U,V,'.k')
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),outXY); hold on
% scatter(U,V,'.k')

vecUV = mean(d.sin(:,:),2);
[vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
vecUV = abs(angle(vecUV));

F = scatteredInterpolant(U(ind),V(ind),vecUV(ind),interp,extrap);
vecXY = F(X,Y);
cMap = redblue(256);
% figure('WindowStyle','docked');
% hIm = imagesc(X(1,:),Y(:,1),vecXY); hold on
% scatter(U,V,'.k')

vecXY(outXY|outUV) = nan;
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),vecXY,hIm.Parent.CLim); hold on
% scatter(U,V,'.k')

sigma = filterSD./mode(diff(X(1,:)));
kernel = fspecial('gauss',2*ceil(2*sigma.*[1 1])+1,sigma);
vecXY = nanconv(vecXY,kernel,'nanout');
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),vecXY,hIm.Parent.CLim); hold on
% scatter(U,V,'.k')


% F = scatteredInterpolant(U(ind),V(ind),vecUV(ind),interp,'linear');
% vecXY = F(X,Y);
% cMap = redblue(256);
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),vecXY,hIm.Parent.CLim); hold on
% scatter(U,V,'.k')

% sigma = filterSD./[mode(diff(X(1,:))) mode(diff(Y(:,1)))];
% if all(sigma)
%     vecXY = imgaussfilt(vecXY,sigma,'FilterDomain','spatial');
% end
% vecXY(outXY|outUV) = nan;
% figure('WindowStyle','docked');
% imagesc(X(1,:),Y(:,1),vecXY,hIm.Parent.CLim); hold on
% scatter(U,V,'.k')

vecXY(outXY|outUV) = level-pi*0.01;

M = contourc(X(1,:),Y(:,1),vecXY,ones(2,1).*level);
cont = cell(0);
contLevel = [];
contN = [];
while ~isempty(M)
    contLevel = [contLevel M(1,1)];
    contN = [contN M(2,1)];
    cont = [cont {M(:,2:1+M(2,1))'}];
    M(:,1:1+M(2,1)) = [];
end
[contN,b] = sort(contN,'descend');
contLevel = contLevel(b);
cont = cont(b);

if p.figOption.verbose>=3 ...
    || (p.figOption.verbose>=2 && p.figOption.subjInd==subjInd)
    visibilityFlag = 'on';
else
    visibilityFlag = 'off';
end

showDelayFovContour(X,Y,vecXY,cMap,U(ind),V(ind),cont,[],visibilityFlag)
% showDelayFovContour(X,Y,vecXYcont,cMap,U(ind),V(ind),cont,[],visibilityFlag)
addEccRef(d,p)
title(['subj' num2str(subjInd)])
drawnow


contData.cont = cont;
contData.contLevel = contLevel;
contData.U = U;
contData.V = V;
contData.vecUV = vecUV;
contData.outUV = outUV;
contData.X = X;
contData.Y = Y;
contData.outXY = outXY;
contData.cMap = cMap;
