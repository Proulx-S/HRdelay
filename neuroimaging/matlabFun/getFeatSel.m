function featSel = getFeatSel(d,p,featSel_fov)
allFeatVal = cell(0);
allFeatP = cell(0);
allFeatMethod = cell(0);
allFeatIndStart = cell(0);
allFeatIndIn = cell(0);

%% Precompute stats on response vector (random effect)
if p.featSel.respVecSig.doIt || p.featSel.respVecDiff.doIt
    statLabel = 'Hotelling'; % 'Pillai' 'Wilks' 'Hotelling' 'Roy'
    condIndPairList = [{[1 2 3]} {[1 2]} {[1 3]} {[2 3]}];
    interceptStat = nan(size(d.sin,1),length(condIndPairList));
    interceptP = nan(size(d.sin,1),length(condIndPairList));
    condStat = nan(size(d.sin,1),length(condIndPairList));
    condP = nan(size(d.sin,1),length(condIndPairList));
    if p.perm.doIt
        parfor condIndPairInd = 1:length(condIndPairList)
            [condStat(:,condIndPairInd),condP(:,condIndPairInd),interceptStat(:,condIndPairInd),interceptP(:,condIndPairInd)] = getDiscrimStats(d,p,condIndPairList{condIndPairInd},statLabel);
        end
    else
        for condIndPairInd = 1:length(condIndPairList)
            [condStat(:,condIndPairInd),condP(:,condIndPairInd),interceptStat(:,condIndPairInd),interceptP(:,condIndPairInd)] = getDiscrimStats(d,p,condIndPairList{condIndPairInd},statLabel);
        end
    end
end



%% Voxels representing stimulus fov
curInfo1 = {'retinoFov'};
featVal = featSel_fov.featVal;
pVal = nan(size(featVal));
% thresh = nan;
curInfo2 = {featSel_fov.featMethod};
startInd = true(size(d.sin,1),1);
curIndIn = logical(featSel_fov.featIndIn);
indFovIn = curIndIn;

allFeatVal(end+1) = {uniformizeOutputs(featVal,condIndPairList)};
allFeatP(end+1) = {uniformizeOutputs(pVal,condIndPairList)};
allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
allFeatIndStart(end+1) = {uniformizeOutputs(startInd,condIndPairList)};
allFeatIndIn(end+1) = {uniformizeOutputs(curIndIn,condIndPairList)};



%% Activated voxels (fixed-effect sinusoidal fit of BOLD timeseries)
% !!!Warning!!! Contrary to what defineFeatSel.m may suggest, activated
% voxels (model explaining the BOLD timeseries) are not detected using only
% the to-be-decoded conditions. It always uses all conditions. The reason
% is that I did not go back to to have processResponses.m <- runGLMs.m
% output stats from fits using only the corresponding conditions. Lazy
% programming, but likely no impact at all.
if p.featSel.act.doIt
    curInfo1 = {'act'};
    
    featVal = d.featSel.F.act;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    startInd = indFovIn;
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                fdr(startInd) = mafdr(pVal(startInd),'BHFDR',true);
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            pVal = featVal.p;
            featVal = featVal.F;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(startInd),curThresh);
        otherwise
            error('X')
    end
    curIndIn = repmat(curIndIn,[1 length(condIndPairList)]);

    allFeatVal(end+1) = {uniformizeOutputs(featVal,condIndPairList)};
    allFeatP(end+1) = {uniformizeOutputs(pVal,condIndPairList)};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndStart(end+1) = {uniformizeOutputs(startInd,condIndPairList)};
    allFeatIndIn(end+1) = {uniformizeOutputs(curIndIn,condIndPairList)};
end


%% Most significant response vectors
if p.featSel.respVecSig.doIt
    curInfo1 = {'respVecSig'};
    featVal = interceptStat; clear interceptStat
    pVal = interceptP; clear interceptP

    featVal;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    startInd = indFovIn;
    switch thresh.threshMethod
        case {'p' 'fdr'}
            pVal;
            featVal;
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                for condIndPairInd = 1:size(pVal,2)
                    fdr(startInd,condIndPairInd) = mafdr(pVal(startInd,condIndPairInd),'BHFDR',true);
                end
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            error('double-check that')
            pVal;
            featVal;
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = featVal>=prctile(featVal(startInd),curThresh);
        otherwise
            error('X')
    end
    allFeatVal(end+1) = {uniformizeOutputs(featVal,condIndPairList)};
    allFeatP(end+1) = {uniformizeOutputs(pVal,condIndPairList)};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndStart(end+1) = {uniformizeOutputs(startInd,condIndPairList)};
    allFeatIndIn(end+1) = {uniformizeOutputs(curIndIn,condIndPairList)};
end


%% Non vein voxels
if p.featSel.vein.doIt
    curInfo1 = {'vein'};
    
    featVal = mean(d.featSel.vein.map(:,:),2);
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    
    % Compute vein %tile threshold on active fov voxels that were identified
    % using all three conditions, since veins should be veins irrespective
    % of which conditions are being decoded)
    startInd = all(catcell(3,allFeatIndIn),3);
    condIndPairList_tmp = condIndPairList;
    for i = 1:length(condIndPairList)
        condIndPairList_tmp{i} = num2str(condIndPairList{i});
    end
    startInd = startInd(:,ismember(condIndPairList_tmp,num2str([1 2 3])));
    switch thresh.threshMethod
        case '%ile'
            pVal = nan(size(featVal));
            featVal;
            curThresh = 100-thresh.percentile;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            curIndIn = featVal<=prctile(featVal(startInd,1),curThresh);            
            
            allFeatVal(end+1) = {uniformizeOutputs(featVal,condIndPairList)};
            allFeatP(end+1) = {uniformizeOutputs(pVal,condIndPairList)};
            allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
            allFeatIndStart(end+1) = {uniformizeOutputs(startInd,condIndPairList)};
            allFeatIndIn(end+1) = {uniformizeOutputs(curIndIn,condIndPairList)};
        otherwise
            error('X')
    end
end

%% Most discrimant voxels (of the activated voxels)
if p.featSel.respVecDiff.doIt
    curInfo1 = {'respVecDiff'};
    featVal = condStat;
    pVal = condP;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    
    switch thresh.threshMethod
        case {'p' 'fdr'}
            error('double-check that')
            curThresh = thresh.threshVal;
            curInfo2 = {[curInfo2{1} '<' num2str(curThresh)]};
            
            if strcmp(thresh.threshMethod,'fdr')
                fdr = nan(size(pVal));
                for condIndPairInd = 1:size(pVal,2)
                    fdr(startInd,condIndPairInd) = mafdr(pVal(startInd,condIndPairInd),'BHFDR',true);
                end
                curIndIn = fdr<=curThresh;
            else
                curIndIn = pVal<=curThresh;
            end
        case '%ile'
            curThresh = thresh.percentile;
            curInfo2 = {[curInfo2{1} '>' num2str(curThresh)]};
            
            curIndIn = false(size(featVal));
            for condIndPairInd = 1:length(condIndPairList)
                % Compute discriminent voxel %tile threshold in a way that
                % follows the use of p.featSel.global.method in
                % runDecoding.m
                [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(allFeatMethod,condIndPairList,p.featSel.global.method,condIndPairInd);
                startInd_nSpec = all(catcell(3,allFeatIndIn(ind_nSpecFeatSel)),3);
                startInd_nSpec = startInd_nSpec(:,ind_nSpecFeatSelCond);
                startInd_spec = all(catcell(3,allFeatIndIn(ind_specFeatSel)),3);
                startInd_spec = startInd_spec(:,ind_specFeatSelCond);
                startInd = startInd_nSpec & startInd_spec;
                
                curIndIn(:,condIndPairInd) = featVal(:,condIndPairInd)>=prctile(featVal(startInd,condIndPairInd),curThresh);
            end
        otherwise
            error('X')
    end
    
    allFeatVal(end+1) = {uniformizeOutputs(featVal,condIndPairList)};
    allFeatP(end+1) = {uniformizeOutputs(pVal,condIndPairList)};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndStart(end+1) = {uniformizeOutputs(startInd,condIndPairList)};
    allFeatIndIn(end+1) = {uniformizeOutputs(curIndIn,condIndPairList)};
end


%% Sumarize feature selection output
for i = 1:size(allFeatVal,2)
    allFeatVal{i} = permute(allFeatVal{i},[1 3 2]);
    allFeatP{i} = permute(allFeatP{i},[1 3 2]);
    allFeatIndStart{i} = permute(allFeatIndStart{i},[1 3 2]);
    allFeatIndIn{i} = permute(allFeatIndIn{i},[1 3 2]);
end

featSel.featSeq.featVal = catcell(2,allFeatVal);
featSel.featSeq.featP = catcell(2,allFeatP);
featSel.featSeq.featQtile = nan(size(featSel.featSeq.featVal));
featSel.featSeq.featIndStart = catcell(2,allFeatIndStart);
featSel.featSeq.featIndIn = catcell(2,allFeatIndIn);
featSel.featSeq.featSelList = allFeatMethod;
featSel.featSeq.condPairList = permute(condIndPairList,[1 3 2]);
featSel.featSeq.info = 'vox X featSel X condPair';
featSel.featSeq.info2 = p.featSel.global.method;

% Compute quantile
for featInd = 1:size(featSel.featSeq.featQtile,2)
    for condIndPairInd = 1:size(featSel.featSeq.featQtile,3)
        x = featSel.featSeq.featVal(:,featInd,condIndPairInd);
        startInd = featSel.featSeq.featIndStart(:,featInd,condIndPairInd);
        %                 startInd = true(size(x));
        [fx,x2] = ecdf(x(startInd));
        x2 = x2(2:end); fx = fx(2:end);
        [~,b] = ismember(x,x2);
        featSel.featSeq.featQtile(b~=0,featInd,condIndPairInd) = fx(b(b~=0));
    end
end

% switch p.featSel.global.method
%     case 'all'
%         indIn = all(featSel.featSeq.featIndIn,2);
%         
%     otherwise
%         error('X')
% end
%% Output
% featSel.indIn = indIn;
% featSel.condPairList = permute(condIndPairList,[1 3 2]);
% featSel.info = 'vox X featSel X condPair';





function [featVal,pVal,featVal2,pVal2] = getDiscrimStats(d,p,condIndPair,statLabel)
if ~exist('statLabel','var') || isempty(statLabel)
    statLabel = 'Hotelling'; % 'Pillai' 'Wilks' 'Hotelling' 'Roy'
end

% Compute stats
voxIndList = 1:size(d.sin,1);
[x,y,~] = getXYK(d,p);
[x,~] = polarSpaceNormalization(x,'cart');

% dummy pass with the manova2.m custom wrapper to obtain design matrix
ind = ismember(y,condIndPair);
withinDesign = table({'real' 'imag'}','VariableNames',{'complex'});
withinModel = 'complex';
voxInd = 1;
t = table(cellstr(num2str(y(ind))),real(x(ind,voxIndList(voxInd))),imag(x(ind,voxIndList(voxInd))),...
    'VariableNames',{'cond','real','imag'});
rm = fitrm(t,'real,imag~cond','WithinDesign',withinDesign,'WithinModel',withinModel);
Xmat = rm.DesignMatrix;
[TBL,A,C_MV,D,withinNames,betweenNames] = manova2(rm,withinModel,[],statLabel);

% actual pass with the manova3.m custom wrapper for fast computation of only
% the relevant stats
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

function featVec = uniformizeOutputs(featVec,condIndPairList)
if size(featVec,2)==1
    featVec = repmat(featVec,[1 length(condIndPairList) 1]);
elseif size(featVec,2)~=length(condIndPairList)
    error('X')
end

