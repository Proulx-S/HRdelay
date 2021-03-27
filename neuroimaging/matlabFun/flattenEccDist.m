function d = flattenEccDist(d,p)
if p.figOption.verbose>1
    plotFlag = 1;
else
    plotFlag = 0;
end

%% Get pdf(ecc) for each subject
sessInd = 1;
for subjInd = 1:size(d,1)
    [eccXY{subjInd},densityXY{subjInd},eccXYspline{subjInd},densityXYspline{subjInd}] = getEccPd(d{subjInd,sessInd},plotFlag);
    if plotFlag
        title(['subj' num2str(subjInd)]);
        figure('WindowStyle','docked');
        hScat = scatter(eccXY{subjInd},densityXY{subjInd},'ko','filled'); hold on
        alpha(hScat,0.1)
        plot(eccXYspline{subjInd},densityXYspline{subjInd},'r')
        title(['subj' num2str(subjInd)]);
    end
end
%% Get the group-averaged pdf(ecc)
ecc = linspace(min([eccXYspline{:}]),max([eccXYspline{:}]),1000);
densityecc = nan(size(d,1),1000);
for subjInd = 1:size(d,1)
    densityecc(subjInd,:) = interp1(eccXYspline{subjInd},densityXYspline{subjInd},ecc);
end
if plotFlag
    figure('WindowStyle','docked');
    hScat = scatter(cat(1,eccXY{:}),cat(1,densityXY{:}),'ko','filled'); hold on
    alpha(hScat,0.03)
    plot(ecc,densityecc,':r')
    plot(ecc,mean(densityecc,1),'-r')
end
densityecc = nanmean(densityecc,1);
ecc = ecc(~isnan(densityecc));
densityecc = densityecc(~isnan(densityecc));

%% Get the group-averaged cdf(ecc)
cdfecc = cumsum(densityecc);
cdfecc = cdfecc./max(cdfecc);

%% Get the ideal cdf(ecc), one that makes pdf(ecc) flat on a circle
cdfecc2 = linspace(0,1,100);
ecc2 = 2*pi*cdfecc2.^2;
ecc2 = ecc2./max(ecc2).*max(ecc);

%% Put each subject on the group-averaged cdf(ecc) then transform ecc to the ideal cdf(ecc)
for subjInd = 1:size(d,1)
    %group-averaged cdf(ecc)
    curecc = d{subjInd,1}.voxProp.ecc;
    curcdfecc = interp1(ecc,cdfecc,curecc);
    %ideal cdf(ecc)
    cureccFlat = interp1(cdfecc2,ecc2,curcdfecc);
    if plotFlag && subjInd == 1
        figure('WindowStyle','docked');
        plot(ecc,cdfecc,'k'); hold on
        plot(curecc,curcdfecc,'ok'); hold on
        plot(ecc2,cdfecc2,'b'); hold on
        plot(cureccFlat,curcdfecc,'ob'); hold on
    end
    d{subjInd,1}.voxProp.eccFlat = cureccFlat;
    d{subjInd,2}.voxProp.eccFlat = cureccFlat;
end
% 
% 
% % ecc2 = [0 max(ecc)*fac];
% % cdfecc2 = [0 1];
% % ecc2 = linspace(min(ecc),max(ecc),length(unique(ecc)));
% % cdfecc2 = linspace(1/length(unique(ecc)),1,length(unique(ecc)));
% % ecc2 = linspace(min(ecc),max(ecc),length(unique(ecc)));
% % cdfecc2 = sqrt(ecc2);
% %% Transform ecc for a flat pdf(ecc) in circle space
% cdfecc2 = linspace(0,1,100);
% ecc2 = 2*pi*cdfecc2.^2;
% ecc2 = ecc2./max(ecc2).*max(ecc);
% if plotFlag
%     plot(ecc2,cdfecc2,'b'); hold on
% end
% for subjInd = 1:size(d,1)
%     vq = interp1(x,v,xq)
%     eccFlat = interp1(cdfecc2,ecc2,cdfecc);
%     d.voxProp.eccFlat = eccFlat;
% end



% if plotFlag
%     figure('WindowStyle','docked');
%     scatter(d{1,1}.voxProp.ecc,d{1,1}.voxProp.eccFlat)
%     xlabel('ecc')
%     ylabel('ecc flat')
% end
% 
% if plotFlag
%     ind = true([size(d{1,1}.sin,1) 1]);
%     plotVoxOnFoV(d{1,1},p,ind,0)
%     title('before')
%     plotVoxOnFoV(d{1,1},p,ind,1)
%     title('after')
% end




function [eccXY,densityXY,eccXYspline,densityXYspline] = getEccPd(d,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 0;
end
ind = true([size(d.sin,1) 1]);
% plotVoxOnFoV(d,p,ind,0)
pol = d.voxProp.pol;
pol(d.voxProp.hemifieldL) = d.voxProp.pol(d.voxProp.hemifieldL)./180*pi;
pol(d.voxProp.hemifieldR) = -d.voxProp.pol(d.voxProp.hemifieldR)./180*pi;
pol = wrapToPi(pol+pi/2);
ecc = d.voxProp.ecc;
[U,V] = pol2cart(pol,ecc);
n = 100;
bwFac = 0.009;
delta = (max([U; V]) - min([U; V]))/n;
[X,Y] = meshgrid(linspace(min([U; V])-delta,max([U; V])+delta,n));
[densityXY,~] = ksdensity([U V],[X(:) Y(:)],'Bandwidth',delta*n*bwFac);
% n = size(xi,1);
% X = xi(1:sqrt(n):end,1);
% Y = xi(1:sqrt(n),2);
% [X,Y] = meshgrid(X,Y);
densityXY = reshape(densityXY,[n n]);

% X = xi(:,1); Y = xi(:,2);
% [X,Y] = meshgrid(xi(:,1),xi(:,2));

% [bandwidth,densityXY,X,Y]=kde2d([U V]);
densityUV = interp2(X,Y,densityXY,U,V);
densityXY(densityXY<(min(densityUV)/10)) = nan;
% densityXY(densityXY<min(densityUV)) = nan;
if plotFlag
    figure('WindowStyle','docked');
    surf(X,Y,densityXY,'LineStyle','none'); hold on
    scatter3(U,V,densityUV,eps,'w.')
    colormap hot
    alpha(.8)
    set(gca, 'color', 'b');
    ax = gca;
    ax.PlotBoxAspectRatio(1:2) = max(ax.PlotBoxAspectRatio(1:2));
    ax.DataAspectRatio(1:2) = max(ax.DataAspectRatio(1:2));
end
[~,eccXY] = cart2pol(X,Y);
eccXY = eccXY(~isnan(densityXY));
densityXY = densityXY(~isnan(densityXY));
sm = 0.8;
[fitresult, ~] = fitSpline(eccXY, densityXY, sm);
eccXYspline = linspace(min(eccXY),max(eccXY),100);
densityXYspline = fitresult(eccXYspline);







% 
% 
% 
% eccXY       = eccXY(~isnan(densityXY));
% densityXY = densityXY(~isnan(densityXY));
% save tmp rho densityXY
% 
% % [tmp,density2] = ecdf(densityXY(:));
% % tmp = [0; diff(tmp)];
% % [~,b] = max(tmp);
% % densityXY(densityXY<(density2(b))*fac) = nan;
% 
% 
% %% compute cdf over sectors
% dSect = 30;
% sectList = 0:dSect:180; sectList(end) = [];
% R = d.voxProp.ecc;
% Rp_hemiL = nan(size(R,1),length(sectList));
% Rp_hemiR = nan(size(R,1),length(sectList));
% for sectInd = 1:length(sectList)
%     indSect = sectList(sectInd)<d.voxProp.pol & d.voxProp.pol<sectList(sectInd)+dSect;
%     
%     curR = d.voxProp.ecc(indSect & d.voxProp.hemifieldL);
%     curRp_hemiL = cdf(nonparamDistFit(curR,0),R);
%     curRp_hemiL(curRp_hemiL==0 | curRp_hemiL==1) = nan;
%     Rp_hemiL(:,sectInd) = curRp_hemiL;
%     
%     curR = d.voxProp.ecc(indSect & d.voxProp.hemifieldR);
%     curRp_hemiR = cdf(nonparamDistFit(curR,0),R);
%     curRp_hemiR(curRp_hemiR==0 | curRp_hemiR==1) = nan;
%     Rp_hemiR(:,sectInd) = curRp_hemiR;
% end
% [~,b] = sort(R);
% imagesc(Rp_hemiL(b,:)')
% 
% curRp_hemiL = curRp_hemiL./length(sectList);
% curRp_hemiR = curRp_hemiR./length(sectList);
% 
% figure('WindowStyle','docked');
% plot(R,curRp_hemiL,'.');
% 
% 
% 
% R = d.voxProp.ecc;
% hist(d.voxProp.pol(d.voxProp.hemifieldL))
% hist(d.voxProp.pol(d.voxProp.hemifieldR))
% Rp = cdf(nonparamDistFit(R,plotFlag),R);
% if plotFlag
%     title('before')
%     xlabel('ecc')
% end
% V = linspace(min(R),max(R)*fac,length(unique(R)));%R
% x = linspace(1/length(unique(R)),1,length(unique(R)));%Rp
% xq = Rp;%Rp
% vq = interp1(x,V,xq);%R
% R2 = vq;
% Rp2 = xq;
% figure('WindowStyle','docked');
% [~,uInd,~] = unique(R);
% plot(R(uInd),Rp(uInd),'.'); hold on
% [~,uInd,~] = unique(R2);
% plot(R2(uInd),Rp2(uInd),'.'); hold on
% nonparamDistFit(R2,plotFlag)
% if plotFlag
%     title('after')
%     xlabel('ecc')
% end
% 
% d.voxProp.eccFlat = R2;
% 
% 
% if plotFlag
%     figure('WindowStyle','docked');
%     scatter(d.voxProp.ecc,d.voxProp.eccFlat)
%     xlabel('ecc')
%     ylabel('ecc flat')
% end
% 
% if plotFlag
%     ind = true([size(d.sin,1) 1]);
%     plotVoxOnFoV(d,p,ind,0)
%     title('before')
%     plotVoxOnFoV(d,p,ind,1)
%     title('after')
% end