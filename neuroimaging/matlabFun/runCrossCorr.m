function [r, rFixDelay] = runCrossCorr(y,p,constantFromRun,plotIt)
if ~isempty(constantFromRun)
    includeConstant = 0;
else
    includeConstant = 1;
end

model = p.model;
fixDelay = p.fixDelay;
t = p.t;
if ~exist('plotIt','var')
    plotIt = 0;
end

% halfModel = 1;
% modelFull = model;
% if halfModel
%     model = model(1:end/2,:);
% end

%Run correlation over all models with no intercept
slope = nan(size(model,1),1);
Rsquared = nan(size(model,1),1);
F = nan(size(model,1),1);
ps = nan(size(model,1),1);
residuals = nan(size(model,1),length(y));

for i = 1:size(model,1)
    if includeConstant
        curModel = [ones(size(model(i,:)')) model(i,:)'];
        cury = y';
        [B,BINT,R,RINT,STATS] = regress(cury,curModel);
        slope(i) = B(2);
        constant(i) = B(1);
        Rsquared(i) = STATS(1);
        F(i) = STATS(2);
        ps(i) = STATS(3);
        residuals(i,:) = R;
    else
        curModel = model(i,:)';
        cury = y'+constantFromRun;
        [B,BINT,R,RINT,STATS] = regress(cury,curModel);
        slope(i) = B(1);
        constant(i) = constantFromRun;
        Rsquared(i) = STATS(1);
        residuals(i,:) = R;
        
% %         [slope(i),~,residuals(i,:),a,b] = regress(y',curModel);
%         %     slope(i) = model(i,:)'\y';
%         Y = slope(i)*model(i,:);
%         RsquaredAdj(i) = getAdjRsquared(y,Y);
%         [F(i),ps(i)] = getF(y,Y,2); %beware, this is wrong !!!
    end
end
residualsNull(1,:) = cury;

% figure('windowStyle','docked');
% plot(t,y); hold on
% plot(t,model'.*repmat(slope',1,120)')
% 
% figure('windowStyle','docked');
% plot(RsquaredAdj)
% figure('windowStyle','docked');
% plot(p)

[~,tmpRmaxInd] = max(Rsquared);
r.lag = tmpRmaxInd;
r.Rsquared = Rsquared(tmpRmaxInd);
r.F = F(tmpRmaxInd);
r.RMSEall = sqrt(nanmean(residuals.^2,2));
r.RMSE = sqrt(nanmean(residuals(tmpRmaxInd,:).^2,2));
r.RMSEnull = sqrt(nanmean(residualsNull.^2,2));
r.p = ps(tmpRmaxInd);
r.slope = slope(tmpRmaxInd);
r.constant = constant(tmpRmaxInd);

[~,tmpRmaxInd] = min(abs(p.tShifts-fixDelay));
rFixDelay.lag = tmpRmaxInd;
rFixDelay.Rsquared = Rsquared(tmpRmaxInd);
rFixDelay.F = F(tmpRmaxInd);
rFixDelay.RMSE = sqrt(nanmean(residuals(tmpRmaxInd,:).^2,2));
rFixDelay.p = ps(tmpRmaxInd);
rFixDelay.slope = slope(tmpRmaxInd);
rFixDelay.constant = constant(tmpRmaxInd);
    
    
%     
%     [tmpRmax,tmpRmaxInd] = max(RsquaredAdj);
%     r.lag = tmpRmaxInd;
%     r.RsquaredAdj = RsquaredAdj(tmpRmaxInd);
%     r.F = F(tmpRmaxInd);
%     r.RMSE = sqrt(nanmean(residuals(tmpRmaxInd,:))^2);
%     r.p = ps(tmpRmaxInd);
%     r.slope = slope(tmpRmaxInd);


% if halfModel
%     %Get maxima at different phase (R)
%     [tmpRmax,tmpRmaxInd] = max(RsquaredAdj(1:size(model,1)/2));
%     r.RsquaredAdj = abs(RsquaredAdj(tmpRmaxInd));
%     r.F = F(tmpRmaxInd);
%     r.p = p(mpRmaxInd);
%     r.slope = slope(tmpRmaxInd);
%     if r.slope<0
%         r.slope = -r.slope;
%         r.lag = tmpRmaxInd+size(model,2)/2;
%     else
%         r.lag = tmpRmaxInd;
%     end
%     
% else
%     
%     %Get maxima at different phase (R)
%     [tmpRmax(1),tmpRmaxInd(1)] = max(RsquaredAdj(1:size(model,1)/2));
%     [tmpRmax(2),tmpRmaxInd(2)] = max(RsquaredAdj(size(model,1)/2+1:end));
%     tmpRmaxInd(2) = tmpRmaxInd(2)+size(model,1)/2;
%     tmpP1 = slope(tmpRmaxInd);
%     [~,b]=max(tmpP1);
%     slectedShift = tmpRmaxInd(b);
%     
%     r.RsquaredAdj = abs(RsquaredAdj(slectedShift));
%     r.F = F(slectedShift);
%     r.p = p(slectedShift);
%     r.lag = slectedShift;
%     r.slope = slope(slectedShift);
%     
% end


% %Get maxima at different phase (abs of R)
% [tmpRmax(1),tmpRmaxInd(1)] = max(abs(RsquaredAdj(1:size(model,1)/2)));
% [tmpRmax(2),tmpRmaxInd(2)] = max(abs(RsquaredAdj(size(model,1)/2+1:end)));
% tmpRmaxInd(2) = tmpRmaxInd(2)+size(model,1)/2;
% tmpP1 = slope(tmpRmaxInd);
% [~,b]=max(tmpP1);
% slectedShift = tmpRmaxInd(b);
% 
% r.RsquaredAdj = abs(RsquaredAdj(slectedShift));
% r.F = F(slectedShift);
% r.p = p(slectedShift);
% r.lag = slectedShift;
% r.slope = slope(slectedShift);






% %Get maxima at different phase
% [~,tmpPminInd(1)] = min(p(1:size(model,1)/2));
% [~,tmpPminInd(2)] = min(p(size(model,1)/2+1:end));
% tmpPminInd(2) = tmpPminInd(2)+size(model,1)/2;
% 
% tmpP1 = slope(tmpPminInd);
% [~,b]=max(tmpP1);
% slectedShift = tmpPminInd(b);
% 
% r.RsquaredAdj = abs(RsquaredAdj(slectedShift));
% r.F = F(slectedShift);
% r.p = p(slectedShift);
% r.lag = slectedShift;
% r.slope = slope(slectedShift);





%Plot
if plotIt
    figure('windowStyle','docked');
    plot(t,y); hold on
    if includeConstant
        plot(t,model(r.lag,:)*r.slope+r.constant,'r')
        plot(t,r.constant,'k')
    else
        plot(t,model(r.lag,:)*r.slope,'r')
    end
    
    figure('windowStyle','docked');
    plot(t,y); hold on
    if includeConstant
        plot(t,model(rFixDelay.lag,:)*rFixDelay.slope+rFixDelay.constant,'r')
        plot(t,rFixDelay.constant,'k')
    else
        plot(t,model(rFixDelay.lag,:)*rFixDelay.slope,'r')
    end
end
