function out = getAllMetric2(rIn)


for anaInd = 1:length(rIn)
    nTrial = length(rIn{anaInd}.yhat(~isnan(rIn{anaInd}.yhat(:,1,1,1)),1,1,1))/2;
    for fold = 1:size(rIn{anaInd}.yhat,3)
        for subFold = 1:size(rIn{anaInd}.yhat,2)
            yhat = squeeze(rIn{anaInd}.yhat(:,subFold,fold,:));
            
            yhat1 = yhat(rIn{anaInd}.y==1);
            yhat1 = reshape(yhat1,[length(find(rIn{anaInd}.y==1))/size(yhat,2) size(yhat,2)]);
            yhat1 = reshape(yhat1(~isnan(yhat1)),nTrial,size(yhat1,2));
            yhat2 = yhat(rIn{anaInd}.y==2);
            yhat2 = reshape(yhat2,[length(find(rIn{anaInd}.y==2))/size(yhat,2) size(yhat,2)]);
            yhat2 = reshape(yhat2(~isnan(yhat2)),nTrial,size(yhat2,2));

            
%             yhat1 = yhat(rIn{anaInd}.y(:,1)==1,:);
%             yhat1 = reshape(yhat1(~isnan(yhat1)),nTrial,size(yhat1,2));
%             yhat2 = yhat(rIn{anaInd}.y(:,1)==2,:);
%             yhat2 = reshape(yhat2(~isnan(yhat2)),nTrial,size(yhat2,2));
            
            
%             yhat = rIn{anaInd}.yhat(rIn{anaInd}.y(:,1)==1,subFold,fold,:);
%             yhat1 = reshape(yhat(~isnan(yhat)),nTrial,size(yhat,4));
%             yhat = rIn{anaInd}.yhat(rIn{anaInd}.y(:,1)==2,subFold,fold,:);
%             yhat2 = reshape(yhat(~isnan(yhat)),nTrial,size(yhat,4));
            %Compute metric
            %acc
            acc(:,subFold,fold) = (sum(yhat1>=0,1) + sum(yhat2<0,1)) ./ (size(yhat1,1)+size(yhat2,1)) * 100;
            %accAbs
            accAbs(:,subFold,fold) = abs(acc(:,subFold,fold)-50)+50;
            %dist
            dist(:,subFold,fold) = mean(yhat1,1)-mean(yhat2,1);
            %distAbs
            distAbs(:,subFold,fold) = abs(dist(:,subFold,fold));
%             %distT
%             [~,~,~,STATS] = ttest2(yhat1,yhat2,'dim',1);
%             distT(:,subFold,fold) = STATS.tstat;
%             %distTabs
%             distTabs(:,subFold,fold) = abs(distT(:,subFold,fold));
            %n
            n(:,subFold,fold) = repmat(size(yhat1,1)+size(yhat2,1),size(yhat1,2),1);
            nCorrect(:,subFold,fold) = acc(:,subFold,fold)/100.*n(:,subFold,fold);
        end
    end
    
    yhat = rIn{anaInd}.yhat;
    y = reshape(rIn{anaInd}.y,[size(rIn{anaInd}.y,1) 1 1 size(rIn{anaInd}.y,2)]);
    y = repmat(y,[1 size(rIn{anaInd}.yhat,2) size(rIn{anaInd}.yhat,3) 1]);
    yhat1 = yhat(y==1);
    yhat1 = reshape(yhat1,[length(find(y==1))/numel(y(1,:,:,:)) size(yhat,2) size(yhat,3) size(yhat,4)]);
    yhat1 = squeeze(nanmean(yhat1,3));
    yhat2 = yhat(y==2);
    yhat2 = reshape(yhat2,[length(find(y==2))/numel(y(1,:,:,:)) size(yhat,2) size(yhat,3) size(yhat,4)]);
    yhat2 = squeeze(nanmean(yhat2,3));
    
%     yhat1 = rIn{anaInd}.yhat(rIn{anaInd}.y(:,1)==1,:,:,:);
%     yhat1 = squeeze(nanmean(yhat1,3));
%     yhat2 = rIn{anaInd}.yhat(rIn{anaInd}.y(:,1)==2,:,:,:);
%     yhat2 = squeeze(nanmean(yhat2,3));
    
    
    %distT
    [~,~,~,STATS] = ttest2(mean(yhat1,2),mean(yhat2,2),'dim',1);
    distT(1,1,1) = STATS.tstat;
    %AUC
    for i = 1:size(y,4)
        [~,~,~,AUC] = perfcurve(y(:,1,1,i),cat(1,yhat1(:,i),yhat2(:,i)),1);
        auc(i,1,1,1) = AUC*100;
    end
    out{anaInd}.yhat1 = yhat1;
    out{anaInd}.yhat2 = yhat2;
    out{anaInd}.auc(1,:) = mean(mean(auc,2),3);
    out{anaInd}.acc(1,:) = mean(mean(acc,2),3);
    out{anaInd}.accAbs(1,:) = mean(mean(accAbs,2),3);
    out{anaInd}.dist(1,:) = mean(mean(dist,2),3);
    out{anaInd}.distAbs(1,:) = mean(mean(distAbs,2),3);
    out{anaInd}.distT(1,:) = mean(mean(distT,2),3);
    out{anaInd}.distTabs(1,:) = nan;
%     out{anaInd}.distTabs(1,:) = mean(mean(distTabs,2),3);
    out{anaInd}.n(1,:) = sum(n(:,1,:),3);
    out{anaInd}.nCorrect(1,:) = sum(nanmean(nCorrect,2),3);
    
    if isfield(rIn{anaInd},'info')
        out{anaInd}.info = rIn{anaInd}.info;
    end
    
    if isfield(rIn{anaInd},'subLevel')
        out{anaInd}.subLevel = getAllMetric2(rIn{anaInd}.subLevel);
    end

end

