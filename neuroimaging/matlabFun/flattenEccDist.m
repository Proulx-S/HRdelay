function flattenEccDist(d,p)
% fac = 1/(2*2*pi);
fac = 10;
if p.figOption.verbose>1
    plotFlag = 1;
else
    plotFlag = 0;
end

ind = true([size(d.sin,1) 1]);
plotVoxOnFoV(d,p,ind,0)
pol = d.voxProp.pol;
pol(d.voxProp.hemifieldL) = d.voxProp.pol(d.voxProp.hemifieldL)./180*pi;
pol(d.voxProp.hemifieldR) = -d.voxProp.pol(d.voxProp.hemifieldR)./180*pi;
pol = wrapToPi(pol+pi/2);
ecc = d.voxProp.ecc;
[u,v] = pol2cart(pol,ecc);
[bandwidth,density,X,Y]=kde2d([u v]);
[tmp,density2] = ecdf(density(:));
tmp = [0; diff(tmp)];
[~,b] = max(tmp);
density(density<(density2(b))*fac) = nan;

figure('WindowStyle','docked');
surf(X,Y,density,'LineStyle','none'); hold on
scatter3(u,v,interp2(X,Y,density,u,v),eps,'w.')
colormap hot
alpha(.8)
set(gca, 'color', 'b');
ax = gca;
ax.PlotBoxAspectRatio(1:2) = max(ax.PlotBoxAspectRatio(1:2));
ax.DataAspectRatio(1:2) = max(ax.DataAspectRatio(1:2));



%% compute cdf over sectors
dSect = 30;
sectList = 0:dSect:180; sectList(end) = [];
R = d.voxProp.ecc;
Rp_hemiL = nan(size(R,1),length(sectList));
Rp_hemiR = nan(size(R,1),length(sectList));
for sectInd = 1:length(sectList)
    indSect = sectList(sectInd)<d.voxProp.pol & d.voxProp.pol<sectList(sectInd)+dSect;
    
    curR = d.voxProp.ecc(indSect & d.voxProp.hemifieldL);
    curRp_hemiL = cdf(nonparamDistFit(curR,0),R);
    curRp_hemiL(curRp_hemiL==0 | curRp_hemiL==1) = nan;
    Rp_hemiL(:,sectInd) = curRp_hemiL;
    
    curR = d.voxProp.ecc(indSect & d.voxProp.hemifieldR);
    curRp_hemiR = cdf(nonparamDistFit(curR,0),R);
    curRp_hemiR(curRp_hemiR==0 | curRp_hemiR==1) = nan;
    Rp_hemiR(:,sectInd) = curRp_hemiR;
end
[~,b] = sort(R);
imagesc(Rp_hemiL(b,:)')

curRp_hemiL = curRp_hemiL./length(sectList);
curRp_hemiR = curRp_hemiR./length(sectList);

figure('WindowStyle','docked');
plot(R,curRp_hemiL,'.');



R = d.voxProp.ecc;
hist(d.voxProp.pol(d.voxProp.hemifieldL))
hist(d.voxProp.pol(d.voxProp.hemifieldR))
Rp = cdf(nonparamDistFit(R,plotFlag),R);
if plotFlag
    title('before')
    xlabel('ecc')
end
v = linspace(min(R),max(R)*fac,length(unique(R)));%R
x = linspace(1/length(unique(R)),1,length(unique(R)));%Rp
xq = Rp;%Rp
vq = interp1(x,v,xq);%R
R2 = vq;
Rp2 = xq;
figure('WindowStyle','docked');
[~,uInd,~] = unique(R);
plot(R(uInd),Rp(uInd),'.'); hold on
[~,uInd,~] = unique(R2);
plot(R2(uInd),Rp2(uInd),'.'); hold on
nonparamDistFit(R2,plotFlag)
if plotFlag
    title('after')
    xlabel('ecc')
end

d.voxProp.eccFlat = R2;


if plotFlag
    figure('WindowStyle','docked');
    scatter(d.voxProp.ecc,d.voxProp.eccFlat)
    xlabel('ecc')
    ylabel('ecc flat')
end

if plotFlag
    ind = true([size(d.sin,1) 1]);
    plotVoxOnFoV(d,p,ind,0)
    title('before')
    plotVoxOnFoV(d,p,ind,1)
    title('after')
end