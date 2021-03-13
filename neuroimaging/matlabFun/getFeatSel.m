function [featSel,featSel2] = getFeatSel(d,p,condPair)
ind = true(size(d.sin,1),1);
info1 = {'V1'};
info2 = {'V1'};
n = nnz(ind);


allFeatVal = cell(0);
allFeatP = cell(0);
allFeatFdrInd = cell(0);
allFeatMethod = cell(0);
% allFeatCondPair = cell(0);
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
            
            fdrInd = false(size(featVal));
            curIndIn = featVal<=prctile(featVal(ind),curThresh);
            ind = ind & curIndIn;
            
            allFeatVal(end+1) = {featVal};
            allFeatP(end+1) = {pVal};
            allFeatFdrInd(end+1) = {fdrInd};
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
                fdrInd = ind;
                fdr(fdrInd) = mafdr(pVal(fdrInd),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                fdrInd = false(size(featVal));
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
            
            fdrInd = false(size(featVal));
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
    allFeatFdrInd(end+1) = {fdrInd};
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
                fdrInd = ind;
                fdr(fdrInd) = mafdr(pVal(fdrInd),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                fdrInd = false(size(featVal));
                curIndIn = pVal<=curThresh;
            end
            ind = ind & curIndIn;

%             allFeatDir(end+1) = {'<'};
%             allFeatPerc(end+1) = {thresh};
        case '%ile'
            pVal;
            featVal;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            fdrInd = false(size(featVal));
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
    allFeatFdrInd(end+1) = {fdrInd};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
    info1(end+1) = curInfo1;
    info2(end+1) = curInfo2;
    n(end+1) = nnz(ind);
end

%% Most discrimant voxels
if p.featSel.discrim.doIt
    curInfo1 = {'discrim'};
    condIndPairList = [{[1 2]} {[1 3]} {[2 3]} {[1 2 3]}];
    curFeatVal = nan(size(ind,1),1,length(condIndPairList));
    curFeatP = nan(size(ind,1),1,length(condIndPairList));
    curFeatFdrInd = false(size(ind,1),1,length(condIndPairList));
%     curFeatMethod = cell(1,1,length(condIndPairList));
%     curFeatCondPair = cell(1,1,length(condIndPairList));
    curFeatIndIn = false(size(ind,1),1,length(condIndPairList));
    for condIndPairInd = 1:length(condIndPairList)
        tic
        [featVal,pVal,featVal2,pVal2] = getDiscrimStats(d,p,condIndPairList{condIndPairInd});
        toc
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
                    fdrInd = ind;
                    fdr(fdrInd) = mafdr(pVal(fdrInd),'BHFDR',true);
                    curIndIn = fdr<=curThresh;
                else
                    fdrInd = false(size(featVal));
                    curIndIn = pVal<=curThresh;
                end
%                 ind = ind & curIndIn;
            case '%ile'
                pVal;
                featVal;
                curThresh = thresh.percentile;
                curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
                
                fdrInd = false(size(featVal));
                curIndIn = featVal>=prctile(featVal(ind),curThresh);
%                 ind = ind & curIndIn;
            otherwise
                error('X')
        end
        
        curFeatVal(:,:,condIndPairInd) = featVal;
        curFeatP(:,:,condIndPairInd) = pVal;
        curFeatFdrInd(:,:,condIndPairInd) = fdrInd;
%         curFeatCondPair(:,:,condIndPairInd) = condIndPairList(condIndPairInd);
        curFeatIndIn(:,:,condIndPairInd) = curIndIn;
    end
    
    
    allFeatVal(end+1) = {curFeatVal};
    allFeatP(end+1) = {curFeatP};
    allFeatFdrInd(end+1) = {curFeatFdrInd};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
%     allFeatCondPair(end+1) = {curFeatCondPair};
    allFeatIndIn(end+1) = {curFeatIndIn};
    
%     info1(end+1) = curInfo1;
%     info2(end+1) = curInfo2;
%     n(end+1) = nnz(ind);
else
    condIndPairList = [];
end

%% Multivariate (combined) feature selection
sz = size(allFeatVal{1},[1 2 3]);
for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)>sz(3)
        sz(3) = size(allFeatVal{i},3);
    end
end
sz(2) = size(allFeatVal,2);

featSel2.featSeq.featVal = nan(sz);
featSel2.featSeq.featP = nan(sz);
featSel2.featSeq.featFdrInd = false(sz);
featSel2.featSeq.featQtile = nan(sz);
featSel2.featSeq.featIndIn = false(sz);
featSel2.featSeq.featSelList = allFeatMethod;
featSel2.featSeq.condPairList = permute(condIndPairList,[1 3 2]);
featSel2.featSeq.info = 'vox X featSel X condPair';
for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)==1
        allFeatVal{i} = repmat(allFeatVal{i},[1 1 sz(3)]);
        allFeatP{i} = repmat(allFeatP{i},[1 1 sz(3)]);
        allFeatFdrInd{i} = repmat(allFeatFdrInd{i},[1 1 sz(3)]);
        allFeatIndIn{i} = repmat(allFeatIndIn{i},[1 1 sz(3)]);
    end 
end
if p.featSel.global.doIt
    featSel2.featSeq.featVal = catcell(2,allFeatVal);
    featSel2.featSeq.featP = catcell(2,allFeatP);
    featSel2.featSeq.featFdrInd = catcell(2,allFeatFdrInd);
    featSel2.featSeq.featQtile = nan(size(featSel2.featSeq.featVal));
    featSel2.featSeq.featIndIn = catcell(2,allFeatIndIn);
    
    switch p.featSel.global.method
        case 'all'
            indIn = all(featSel2.featSeq.featIndIn,2);
            % Compute quantile
            for featInd = 1:sz(2)
                for condIndPairInd = 1:sz(3)
                    x = featSel2.featSeq.featVal(:,featInd,condIndPairInd);
                    [fx,x2] = ecdf(x);
                    x2 = x2(2:end); fx = fx(2:end);
                    [~,b] = ismember(x,x2);
                    featSel2.featSeq.featQtile(:,featInd,condIndPairInd) = fx(b);
                end
            end
            
            
            %% Output2
            featSel2.indIn = indIn;
            featSel2.info = 'vox X featSel X condPair';
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


function [featVal,pVal,featVal2,pVal2] = getDiscrimStats(d,p,condIndPair)
% Compute stats
featVal = nan(size(d.sin,1),1);
pVal = nan(size(d.sin,1),1);
featVal2 = nan(size(d.sin,1),1);
pVal2 = nan(size(d.sin,1),1);
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
[TBL,A,C_MV,D,withinNames,betweenNames] = manova2(rm,withinModel);
for voxInd = 1:length(voxIndList)
    Ymat = [real(x(ind,voxIndList(voxInd))) imag(x(ind,voxIndList(voxInd)))];
    [stat,PVAL] = manova3(Xmat,Ymat,C_MV,A,D);
    
    featVal(voxIndList(voxInd)) = stat(ismember(betweenNames,'cond'));
    pVal(voxIndList(voxInd)) = PVAL(ismember(betweenNames,'cond'));
    
    featVal2(voxIndList(voxInd)) = stat(ismember(betweenNames,'(Intercept)'));
    pVal2(voxIndList(voxInd)) = PVAL(ismember(betweenNames,'(Intercept)'));
end

function [Value,pValue,ds] = getStats(X,A,B,C,D,SSE,withinNames,betweenNames)

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

Value = [L];
F = [L_F];
df1 = [L_df1];
df2 = [L_df2];
pValue = fcdf(F, df1, df2, 'upper');

if exist('withinNames','var') && ~isempty(withinNames) && exist('betweenNames','var') && ~isempty(betweenName)
    Within = withinNames;
    Between = betweenNames;
    Statistic = categorical({'Wilks'}');
    Value = [L];
    F = [L_F];
    RSquare = [L_rsq];
    df1 = [L_df1];
    df2 = [L_df2];
    ds = table(Within,Between,Statistic, Value, F, RSquare,df1,df2);
    ds.pValue = fcdf(F, df1, df2, 'upper');
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

function [stat,PVAL] = manova3(Xmat,Ymat,C,A,D)
[Beta,DFE,Cov] = fitrm2(Xmat,Ymat);
SSE = DFE * Cov;

stat = nan(size(A,1),1);
PVAL = nan(size(A,1),1);
i = 0;
for withinTestInd = 1:length(C)
    for betweenTestInd = 1:length(A)
        i = i+1;
        [stat(i),PVAL(i),~] = getStats(Xmat,A{betweenTestInd},Beta,C{withinTestInd},D,SSE);
    end
end
