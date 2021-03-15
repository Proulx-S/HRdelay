function featSel = getFeatSel(d,p)
ind = true(size(d.sin,1),1);

allFeatVal = cell(0);
allFeatP = cell(0);
allFeatPrevInd = cell(0);
allFeatMethod = cell(0);
allFeatIndIn = cell(0);

%% Non vein voxels
if p.featSel.vein.doIt
    curInfo1 = {'vein'};
    
    featVal = mean(d.featSel.vein.map(:,:),2);
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    prevInd = true(size(featVal));
    switch thresh.threshMethod
        case '%ile'
            pVal = nan(size(featVal));
            featVal;
            curThresh = 100-thresh.percentile;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            curIndIn = featVal<=prctile(featVal(ind),curThresh);
            
            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatPrevInd(end+1) = {prevInd};
            allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
            allFeatIndIn(end+1) = {curIndIn};
        otherwise
            error('X')
    end
    ind = ind & curIndIn;
    
%     info1(end+1) = curInfo1;
%     info2(end+1) = curInfo2;
%     n(end+1) = nnz(ind);
end

%% Activated voxels (fixed-effect sinusoidal fit of BOLD timeseries)
if p.featSel.act.doIt
    curInfo1 = {'act'};
    
    featVal = d.featSel.F.act;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    prevInd = ind;
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                fdr(prevInd) = mafdr(pVal(prevInd),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(prevInd),curThresh);
        otherwise
            error('X')
    end
    ind = ind & curIndIn;

    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatPrevInd(end+1) = {prevInd};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};

%     info1(end+1) = curInfo1;
%     info2(end+1) = curInfo2;
%     n(end+1) = nnz(ind);
end


%% Stats on response vector (random effect)
if p.featSel.respVecSig.doIt || p.featSel.respVecDiff.doIt
    statLabel = 'Hotelling'; % 'Pillai' 'Wilks' 'Hotelling' 'Roy'
    condIndPairList = [{[1 2]} {[1 3]} {[2 3]} {[1 2 3]}];
    interceptStat = nan(size(d.sin,1),length(condIndPairList));
    interceptP = nan(size(d.sin,1),length(condIndPairList));
    condStat = nan(size(d.sin,1),length(condIndPairList));
    condP = nan(size(d.sin,1),length(condIndPairList));
    for condIndPairInd = 1:length(condIndPairList)
        [condStat(:,condIndPairInd),condP(:,condIndPairInd),interceptStat(:,condIndPairInd),interceptP(:,condIndPairInd)] = getDiscrimStats(d,p,condIndPairList{condIndPairInd},statLabel);
    end
else
    condIndPairList = [];
end

%% Most significant response vectors
if p.featSel.respVecSig.doIt
    curInfo1 = {'respVecSig'};
    featVal = interceptStat; clear interceptStat
    pVal = interceptP; clear interceptP

    featVal;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    prevInd = ind;
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal;
            featVal;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                for condIndPairInd = 1:size(pVal,2)
                    fdr(prevInd,condIndPairInd) = mafdr(pVal(prevInd,condIndPairInd),'BHFDR',true);
                end
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            pVal;
            featVal;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(prevInd),curThresh);
        otherwise
            error('X')
    end
    ind = ind & curIndIn;
    
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatPrevInd(end+1) = {prevInd};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
%     info1(end+1) = curInfo1;
%     info2(end+1) = curInfo2;
%     n(end+1) = nnz(ind);
end



% %% Response vector distribution across voxels
% [bandwidth,density,X,Y]=kde2d([real(x') imag(x')]);
% figure('WindowStyle','docked');
% imagesc(X(1,1:end),Y(1:end,1),density);
% grid on
% ax = gca;
% ax.YDir = 'normal';
% figure('WindowStyle','docked');
% scatter(real(x'),imag(x'),'k.')
% ax2 = gca;
% ax2.XLim = ax.XLim;
% ax2.YLim = ax.YLim;
% grid on
% 
% 
% clear all
% X=rand(1000,1); Y=sin(X*10*pi)+randn(size(X))/3; data=[X,Y];
% % apply routine
% [bandwidth,density,X,Y]=kde2d(data);
% % plot the data and the density estimate
% surf(X,Y,density,'LineStyle','none'), view([0,70])
% colormap hot, hold on, alpha(.8)
% set(gca, 'color', 'blue');
% plot(data(:,1),data(:,2),'w.','MarkerSize',5)





%% Most discrimant voxels
if p.featSel.respVecDiff.doIt
    curInfo1 = {'respVecDiff'};
    % Threshold
    featVal = condStat;
    pVal = condP;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    prevInd = ind;                
    switch thresh.threshMethod
        case {'p' 'fdr'}
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                for condIndPairInd = 1:size(pVal,2)
                    fdr(prevInd(:,condIndPairInd),condIndPairInd) = mafdr(pVal(prevInd(:,condIndPairInd),condIndPairInd),'BHFDR',true);
                end
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(prevInd),curThresh);
        otherwise
            error('X')
    end
    ind = ind & curIndIn;
    
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatPrevInd(end+1) = {prevInd};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
%     info1(end+1) = curInfo1;
%     info2(end+1) = curInfo2;
%     n(end+1) = nnz(ind);
end


%% Multivariate (combined) feature selection
for i = 1:size(allFeatVal,2)
    allFeatVal{i} = permute(allFeatVal{i},[1 3 2]);
    allFeatP{i} = permute(allFeatP{i},[1 3 2]);
    allFeatPrevInd{i} = permute(allFeatPrevInd{i},[1 3 2]);
    allFeatIndIn{i} = permute(allFeatIndIn{i},[1 3 2]);
end
sz = size(allFeatVal{1},[1 2 3]);
for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)>sz(3)
        sz(3) = size(allFeatVal{i},3);
    end
end
sz(2) = size(allFeatVal,2);

for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)==1
        allFeatVal{i} = repmat(allFeatVal{i},[1 1 sz(3)]);
    end
    if size(allFeatP{i},3)==1
        allFeatP{i} = repmat(allFeatP{i},[1 1 sz(3)]);
    end
    if size(allFeatPrevInd{i},3)==1
        allFeatPrevInd{i} = repmat(allFeatPrevInd{i},[1 1 sz(3)]);
    end
    if size(allFeatIndIn{i},3)==1
        allFeatIndIn{i} = repmat(allFeatIndIn{i},[1 1 sz(3)]);
    end
end
featSel.featSeq.featVal = catcell(2,allFeatVal);
featSel.featSeq.featP = catcell(2,allFeatP);
featSel.featSeq.featPrevInd = catcell(2,allFeatPrevInd);
featSel.featSeq.featQtile = nan(size(featSel.featSeq.featVal));
featSel.featSeq.featIndIn = catcell(2,allFeatIndIn);
featSel.featSeq.featSelList = allFeatMethod;
featSel.featSeq.condPairList = permute(condIndPairList,[1 3 2]);
featSel.featSeq.info = 'vox X featSel X condPair';


switch p.featSel.global.method
    case 'all'
        indIn = all(featSel.featSeq.featIndIn,2);
        % Compute quantile
        for featInd = 1:sz(2)
            for condIndPairInd = 1:sz(3)
                x = featSel.featSeq.featVal(:,featInd,condIndPairInd);
                [fx,x2] = ecdf(x);
                x2 = x2(2:end); fx = fx(2:end);
                [~,b] = ismember(x,x2);
                featSel.featSeq.featQtile(:,featInd,condIndPairInd) = fx(b);
            end
        end
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
%% Output
featSel.indIn = indIn;
featSel.info = 'vox X featSel X condPair';





function [featVal,pVal,featVal2,pVal2] = getDiscrimStats(d,p,condIndPair,statLabel)
if ~exist('statLabel','var') || isempty(statLabel)
    statLabel = 'Hotelling'; % 'Pillai' 'Wilks' 'Hotelling' 'Roy'
end

% Compute stats
%     voxIndList = find(ind);
voxIndList = 1:size(d.sin,1);
[x,y,~] = getXYK(d,p);
%             [x,~] = polarSpaceNormalization(x,p.svmSpace);

%             [~,~,~,STATS] = ttest(real(x(y==1,ind)),real(x(y==2,ind)));
%             featVal(ind) = abs(STATS.tstat);
% Get stats for H0: no difference between conditions
ind = ismember(y,condIndPair);
withinDesign = table({'real' 'imag'}','VariableNames',{'complex'});
withinModel = 'complex';
voxInd = 1;
t = table(cellstr(num2str(y(ind))),real(x(ind,voxIndList(voxInd))),imag(x(ind,voxIndList(voxInd))),...
    'VariableNames',{'cond','real','imag'});
rm = fitrm(t,'real,imag~cond','WithinDesign',withinDesign,'WithinModel',withinModel);
Xmat = rm.DesignMatrix;
[TBL,A,C_MV,D,withinNames,betweenNames] = manova2(rm,withinModel,[],statLabel);

featVal = nan(size(d.sin,1),1);
pVal = nan(size(d.sin,1),1);
featVal2 = nan(size(d.sin,1),1);
pVal2 = nan(size(d.sin,1),1);

Ymat = permute(x(ind,:),[1 3 2]);
Ymat = cat(2,real(Ymat),imag(Ymat));
for voxInd = voxIndList
    [stat,PVAL] = manova3(Xmat,Ymat(:,:,voxInd),C_MV,A,D,statLabel);
    
    featVal(voxInd) = stat(ismember(betweenNames,'cond'));
    pVal(voxInd) = PVAL(ismember(betweenNames,'cond'));
    
    featVal2(voxInd) = stat(ismember(betweenNames,'(Intercept)'));
    pVal2(voxInd) = PVAL(ismember(betweenNames,'(Intercept)'));
end

function [Value,pValue,ds] = getStats(X,A,B,C,D,SSE,statLabel,withinNames,betweenNames)
% Adapted from RepeatedMeasuresModel.m

% Hypothesis matrix H
% H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
% q = rank(Z);
[H,q] = makeH(A,B,C,D,X);

% Error matrix E
E = C'*SSE*C;

p = rank(E+H);
s = min(p,q);
v = size(X,1)-rank(X);
if p^2+q^2>5
    t = sqrt( (p^2*q^2-4) / (p^2+q^2-5));
else
    t = 1;
end
u = (p*q-2)/4;
r = v - (p-q+1)/2;
m = (abs(p-q)-1)/2;
n = (v-p-1)/2;

switch statLabel
    case 'Wilks'
        % ~~~ Wilks' Lambda = L
        % Formally, L = |E| / |H+E|, but it is more convenient to compute it using
        % the eigenvalues from a generalized eigenvalue problem
        lam = eig(H,E);
        mask = (lam<0) & (lam>-100*eps(max(abs(lam))));
        lam(mask) = 0;
        L_df1 = p*q;
        L_df2 = r*t-2*u;
        if isreal(lam) && all(lam>=0) && L_df2>0
            L = prod(1./(1+lam));
        else
            L = NaN;
            L_df2 = max(0,L_df2);
        end
        L1 = L^(1/t);
        L_F = ((1-L1) / L1) * (r*t-2*u)/(p*q);
        L_rsq = 1-L1;
        
        Value = L;
        F = L_F;
        RSquare = L_rsq;
        df1 = L_df1;
        df2 = L_df2;
    case 'Pillai'
            error('code that')
    case 'Hotelling'
        lam = eig(H,E);
        mask = (lam<0) & (lam>-100*eps(max(abs(lam))));
        lam(mask) = 0;
        
        % ~~~ Hotelling-Lawley trace = U
        % U = trace( H * E^-1 ) but it can also be written as the sum of the
        % eigenvalues that we already obtained above
        if isreal(lam) && all(lam>=0)
            U = sum(lam);
        else
            U = NaN;
            n(n<0) = NaN;
        end
        b = (p+2*n)*(q+2*n) / (2*(2*n+1)*(n-1));
        c = (2 + (p*q+2)/(b-1))/(2*n);
        if n>0
            U_F = (U/c) * (4+(p*q+2)/(b-1)) / (p*q);
        else
            U_F = U * 2 * (s*n+1) / (s^2 * (2*m+s+1));
        end
        U_rsq = U / (U+s);
        U_df1 = s*(2*m+s+1);
        U_df2 = 2*(s*n+1);
        
        Value = U;
        F = U_F;
        RSquare = U_rsq;
        df1 = U_df1;
        df2 = U_df2;
    case 'Roy'
            error('code that')
    otherwise
        error('x')
end
pValue = fcdf(F, df1, df2, 'upper');
if exist('withinNames','var') && ~isempty(withinNames) && exist('betweenNames','var') && ~isempty(betweenName)
    Within = withinNames;
    Between = betweenNames;
    Statistic = categorical({statLabel}');
    ds = table(Within,Between,Statistic, Value, F, RSquare,df1,df2);
    ds.pValue = pValue;
else
    ds = [];
end



function [H,q] = makeH(A,B,C,D,X) % Make hypothesis matrix H
% H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
d = A*B*C - D;
[~,RX] = qr(X,0);
XA = A/RX;
Z = XA*XA';
H = d'*(Z\d);  % note Z is often scalar or at least well-conditioned

if nargout>=2
    q = rank(Z);
end

function [Bmat,dfe,Covar] = fitrm2(Xmat,Ymat)
% Use the design matrix to carry out the fit, and compute
% information needed to store in the object
opt.RECT = true;
Bmat = linsolve(Xmat,Ymat,opt);

Resid = Ymat-Xmat*Bmat;
dfe = size(Xmat,1)-size(Xmat,2);
if dfe>0
    Covar = (Resid'*Resid)/dfe;
else
    Covar = NaN(size(Resid,2));
    
    % Diagnose the issue
    if isscalar(formula.PredictorNames) && coefTerms(1)==1 && all(coefTerms(2:end)==2)
        % It appears there is a single predictor that is
        % basically a row label, so try to offer a helpful
        % message
        warning(message('stats:fitrm:RowLabelPredictor',formula.PredictorNames{1}));
    else
        % Generic message
        warning(message('stats:fitrm:NoDF'));
    end
end

function [stat,PVAL] = manova3(Xmat,Ymat,C,A,D,statLabel)
[Beta,DFE,Cov] = fitrm2(Xmat,Ymat);
SSE = DFE * Cov;

stat = nan(size(A,1),1);
PVAL = nan(size(A,1),1);
i = 0;
for withinTestInd = 1:length(C)
    for betweenTestInd = 1:length(A)
        i = i+1;
        [stat(i),PVAL(i),~] = getStats(Xmat,A{betweenTestInd},Beta,C{withinTestInd},D,SSE,statLabel);
    end
end