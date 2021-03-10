function [featSel,featSel2] = getFeatSel(d,p,subjInd,sessInd)
ind = true(size(d.sin,1),1);
info1 = {'V1'};
info2 = {'V1'};
n = nnz(ind);

allFeatVal = cell(0);
allFeatP = cell(0);
allFeatDir = cell(0);
allFeatPerc = cell(0);

%% Non vein voxels
if p.featSel.vein.doIt
    featVal = mean(d.featSel.vein.map(:,:),2);
    thresh = p.featSel.vein;
    curInfo1 = {'veinScore'};
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case '%ile'
            pVal = nan(size(featVal));
            featVal;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '<' num2str(thresh)]};
            
            ind = ind & featVal<=prctile(featVal(ind),thresh);
            
            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatDir(end+1) = {'<'};
            allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Activated voxels (fixed-effect sinusoidal fit)
if p.featSel.act.doIt
    featVal = d.featSel.F.act;
    thresh = p.featSel.act;
    curInfo1 = {'act'};
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case 'fdr'
            error('double-check that')
            P = featVal.p;
            FDR = nan(size(P)); FDR(ind) = mafdr(P(ind),'BHFDR',true);
            
            featVal = FDR;
            thresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} num2str(thresh)]};
            
            ind = ind & featVal<=thresh;
        case '%ile'
            pVal = featVal.p;
            featVal = featVal.F;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(thresh)]};
            
            ind = ind & featVal>=prctile(featVal(ind),thresh);

            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatDir(end+1) = {'>'};
            allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Most significant response vectors
if p.featSel.respVectSig.doIt
    % Compute stats
    pTmp = p;
    pTmp.condPair = 'all';
    [x,y,~] = getXYK(d,pTmp); clear pTmp
    %remove condition differences
    yList = unique(y);
    xMean = mean(x,1);
    for i = 1:length(yList)
        x(y==yList(i),:) = x(y==yList(i),:) - mean(x(y==yList(i),:),1);
    end
    x = x + xMean; clear xMean
    %get stats for H0: vector length=0
    featVal = nan(size(ind));
    P = nan(size(ind));
%     voxIndList = find(ind);
    voxIndList = 1:length(ind);
    for voxInd = 1:length(voxIndList)
        tmp = cat(2,real(x(:,voxIndList(voxInd))),imag(x(:,voxIndList(voxInd))));
        [~,~,featVal(voxIndList(voxInd)),~,~,~,P(voxIndList(voxInd))] = T2Hot1(tmp);
    end
    thresh = p.featSel.respVectSig;
    curInfo1 = {'sigVec'};
    curInfo2 = {thresh.threshMethod};
    % Threshold
    switch thresh.threshMethod
        case 'fdr'
            error('double-check that')
            FDR = nan(size(P)); FDR(ind) = mafdr(P(ind),'BHFDR',true);
            featVal = FDR;
            thresh = p.featSel.respVectSig.threshVal;
            curInfo1 = {'vecSig'};
            curInfo2 = {['FDR<' num2str(thresh)]};
            ind = ind & featVal<=thresh;
        case '%ile'
            pVal = P;
            featVal;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(thresh)]};
            
            ind = ind & featVal>=prctile(featVal(ind),thresh);

            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatDir(end+1) = {'>'};
            allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Most discrimant voxels
if p.featSel.discrim.doIt
    % Compute stats
    featVal = nan(size(ind));
    P = nan(size(ind));
%     voxIndList = find(ind);
    voxIndList = 1:length(ind);
    [x,y,k] = getXYK(d,p);
    %             [x,~] = polarSpaceNormalization(x,p.svmSpace);
    
    %             [~,~,~,STATS] = ttest(real(x(y==1,ind)),real(x(y==2,ind)));
    %             featVal(ind) = abs(STATS.tstat);
    % Get stats for H0: no difference between conditions
    for voxInd = 1:length(voxIndList)
        tmp = cat(1,...
            cat(2,real(x(y==1,voxIndList(voxInd))),imag(x(y==1,voxIndList(voxInd)))),...
            cat(2,real(x(y==2,voxIndList(voxInd))),imag(x(y==2,voxIndList(voxInd))))...
            );
        stats = T2Hot2d(tmp);
        featVal(voxIndList(voxInd)) = stats.T2;
        P(voxIndList(voxInd)) = stats.P;
    end
    % Threshold
    thresh = p.featSel.discrim;
    curInfo1 = {'discrim'};
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case '%ile'
            pVal = P;
            featVal;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(thresh)]};
            
            ind = ind & featVal>=prctile(featVal(ind),thresh);

            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatDir(end+1) = {'>'};
            allFeatPerc(end+1) = {thresh};
        case 'fdr'
            error('double-check that')
            FDR = nan(size(P)); FDR(ind) = mafdr(P(ind),'BHFDR',true);
            featVal = FDR;
            thresh = p.featSel.discrim.threshVal;
            curInfo1 = {'discrim'};
            curInfo2 = {['FDR<' num2str(thresh)]};
            ind = ind & featVal<=thresh;
        otherwise
            error('X')
    end
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Multivariate feature selection
if p.featSel.global.doIt
    allFeatVal = catcell(2,allFeatVal);
    x = zscore(-log(allFeatVal(:,1)));
    y = zscore(log(allFeatVal(:,2)));
    z = zscore(log(allFeatVal(:,3)));
    
    
    [fx,x2] = ecdf(x);
    x2 = x2(2:end); fx = fx(2:end);
    [~,b] = ismember(x,x2);
    fx = fx(b); clear x2
    
    [fy,y2] = ecdf(y);
    y2 = y2(2:end); fy = fy(2:end);
    [~,b] = ismember(y,y2);
    fy = fy(b); clear y2
    
    [fz,z2] = ecdf(z);
    z2 = z2(2:end); fz = fz(2:end);
    [~,b] = ismember(z,z2);
    fz = fz(b); clear z2
    
%     fyz = fy.*fz;
%     fyz = fyz - min(fyz);
%     fyz = fyz./max(fyz);
    fxyz = fx.*fy.*fz;
    fxyz = fxyz - min(fxyz);
    fxyz = fxyz./max(fxyz);
    
    perc = p.featSel.global.percentile;
%     ind2 = fyz>prctile(fyz,perc) & fx>prctile(fx,perc);
    ind2 = fxyz>prctile(fxyz,perc) & fx>prctile(fx,perc);
    
%     if p.figOption.verbose>=1 && subjInd==p.figOption.subjInd && sessInd==p.figOption.sessInd
%         figure('WindowStyle','docked');
%         scatter3(x(ind2),y(ind2),z(ind2),'k.'); hold on
%         scatter3(x(~ind2),y(~ind2),z(~ind2),'r.');
%         % ax = gca;
%         % ax.CameraPosition = camPos;
%         xlabel(['zscore(-log(' info1{1+1} '))'])
%         ylabel(['zscore(log(' info1{1+2} '))'])
%         zlabel(['zscore(log(' info1{1+3} '))'])
%         legend({'include' 'exclude'})
%         % ax = gca; camPos = ax.CameraPosition;
%         % ax = gca; ax.CameraPosition = camPos;
%     end
    
    ind = ind2;
    n = [length(ind) nnz(ind)];
    n = round(n(2)/n(1)*100);
    info2 = {''};
%     if p.figOption.verbose>=1 && subjInd==p.figOption.subjInd && sessInd==p.figOption.sessInd
%         title(['perc=' num2str(p.featSel.global.percentile) '; ' num2str(n) '% of V1 included'])
%         drawnow
%     end
end

%% Output
featSel.ind = ind;
featSel.info1 = info1;
featSel.info2 = info2;
featSel.info3 = cellstr(num2str(n'));
featSel.n = n;

%% Output2
featSel2.indIn = ind2;
featSel2.perc = fxyz;
featSel2.featVal = allFeatVal;
featSel2.featP = allFeatP;
featSel2.featLabel = info1(2:end);
featSel2.info = 'vox X featSel';

% featSel.ind = ind;
% featSel.info1 = info1;
% featSel.info2 = info2;
% featSel.info3 = cellstr(num2str(n'));
% featSel.n = n;
