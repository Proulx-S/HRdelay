function [featSel,featSel2] = getFeatSel(d,p,subjInd,sessInd)
ind = true(size(d.sin,1),1);
info1 = {'V1'};
info2 = {'V1'};
n = nnz(ind);


allFeatVal = cell(0);
allFeatP = cell(0);
allFeatMethod = cell(0);
allFeatIndIn = cell(0);
% allFeatVal = cell(0);
% allFeatP = cell(0);
% allFeatDir = cell(0);
% allFeatPerc = cell(0);

%% Non vein voxels
if p.featSel.vein.doIt
    curInfo1 = {'vein'};
    
    featVal = mean(d.featSel.vein.map(:,:),2);
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case '%ile'
            pVal = nan(size(featVal));
            featVal;
            curThresh = 100-thresh.percentile;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            curIndIn = featVal<=prctile(featVal(ind),curThresh);
            ind = ind & curIndIn;
            
            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
            allFeatIndIn(end+1) = {curIndIn};
            
%             allFeatVal(end+1) = {featVal};
%             allFeatP(end+1) = {pVal};
%             allFeatDir(end+1) = {'<'};
%             allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Activated voxels (fixed-effect sinusoidal fit)
if p.featSel.act.doIt
    curInfo1 = {'act'};
    
    featVal = d.featSel.F.act;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                fdr(ind) = mafdr(pVal(ind),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
            ind = ind & curIndIn;
            
%             allFeatDir(end+1) = {'<'};
%             allFeatPerc(end+1) = {thresh};
        case '%ile'
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(ind),curThresh);
            ind = ind & curIndIn;

%             allFeatVal(end+1) = {featVal};
%             allFeatP(end+1) = {pVal};
%             allFeatMethod(end+1) = curInfo2;
% %             allFeatDir(end+1) = {'>'};
% %             allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};

    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Most significant response vectors
if p.featSel.respVecSig.doIt
    curInfo1 = {'respVecSig'};
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
    pVal = nan(size(ind));
%     voxIndList = find(ind);
    voxIndList = 1:length(ind);
    for voxInd = 1:length(voxIndList)
        tmp = cat(2,real(x(:,voxIndList(voxInd))),imag(x(:,voxIndList(voxInd))));
        [~,~,featVal(voxIndList(voxInd)),~,~,~,pVal(voxIndList(voxInd))] = T2Hot1(tmp);
    end
    
    % Threshold
    featVal;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal;
            featVal;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                fdr(ind) = mafdr(pVal(ind),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
            ind = ind & curIndIn;

%             allFeatDir(end+1) = {'<'};
%             allFeatPerc(end+1) = {thresh};
        case '%ile'
            pVal;
            featVal;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(ind),curThresh);
            ind = ind & curIndIn;

%             allFeatVal(end+1) = {featVal};
%             allFeatP(end+1) = {pVal};
%             allFeatMethod(end+1) = curInfo2;
%             allFeatIndIn(end+1) = {curIndIn};
% %             allFeatDir(end+1) = {'>'};
% %             allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Most discrimant voxels
if p.featSel.discrim.doIt
    curInfo1 = {'discrim'};
    % Compute stats
    featVal = nan(size(ind));
    pVal = nan(size(ind));
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
        pVal(voxIndList(voxInd)) = stats.P;
    end
    
    % Threshold
    featVal;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal;
            featVal;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                fdr(ind) = mafdr(pVal(ind),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
            ind = ind & curIndIn;

%             allFeatDir(end+1) = {'<'};
%             allFeatPerc(end+1) = {thresh};
        case '%ile'
            pVal;
            featVal;
            thresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(ind),curThresh);
            ind = ind & curIndIn;

%             allFeatVal(end+1) = {featVal};
%             allFeatP(end+1) = {pVal};
%             allFeatMethod(end+1) = curInfo2;
%             allFeatIndIn(end+1) = {curIndIn};
% %             allFeatDir(end+1) = {'>'};
% %             allFeatPerc(end+1) = {thresh};
        otherwise
            error('X')
    end
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Multivariate (combined) feature selection
if p.featSel.global.doIt
    featSel2.featSeq.featVal = catcell(2,allFeatVal);
    featSel2.featSeq.featP = catcell(2,allFeatP);
    featSel2.featSeq.featInfo = allFeatMethod;
    featSel2.featSeq.featIndIn = catcell(2,allFeatIndIn);
    featSel2.featSeq.info = 'vox X featSel';
    
    switch p.featSel.global.method
        case 'all'
            indIn = all(featSel2.featSeq.featIndIn,2);
            
            %% Output2
            featSel2.indIn = indIn;
            featSel2.info = 'vox X featSel';
            featSel2.featSeq = featSel2;
        otherwise
            error('code that')
            
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
    end
    
    ind;
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


% featSel.ind = ind;
% featSel.info1 = info1;
% featSel.info2 = info2;
% featSel.info3 = cellstr(num2str(n'));
% featSel.n = n;
