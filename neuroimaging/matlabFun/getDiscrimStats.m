function [featVal,pVal,featVal2,pVal2] = getDiscrimStats(d,p,condIndPair,statLabel)
if ~exist('statLabel','var') || isempty(statLabel)
    statLabel = 'Hotelling'; % 'Pillai' 'Wilks' 'Hotelling' 'Roy'
end

% Compute stats
if isstruct(d)
    voxIndList = 1:size(d.sin,1);
    [x,y,~] = getXYK(d,p);
    [x,~] = polarSpaceNormalization(x,'cart');
    
    featVal = nan(size(d.sin,1),1);
    pVal = nan(size(d.sin,1),1);
    featVal2 = nan(size(d.sin,1),1);
    pVal2 = nan(size(d.sin,1),1);

else
    if size(d,3)>2
        error('X')
    end
    x = d;
    condIndPair = 1:size(x,1);
    voxIndList = 1;
    
    % Remove random-effect of subject
    rho = abs(x) ./ abs(mean(x,1)) .* abs(mean(x(:)));
    theta = wrapToPi(angle(x) - angle(mean(x,1)) + angle(mean(x(:))));
    [u,v] = pol2cart(theta,rho);
    x = complex(u,v);

    % Average sessions
    x = mean(x,3);
    
    y = condIndPair;
    y = repmat(y,[size(x,2) 1])';
    y = y(:);
    x = x(:);
end

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

Ymat = permute(x(ind,:),[1 3 2]);
Ymat = cat(2,real(Ymat),imag(Ymat));
for voxInd = voxIndList
    [stat,PVAL] = manova3(Xmat,Ymat(:,:,voxInd),C_MV,A,D,statLabel);
    
    featVal(voxInd) = stat(ismember(betweenNames,'cond'));
    pVal(voxInd) = PVAL(ismember(betweenNames,'cond'));
    
    featVal2(voxInd) = stat(ismember(betweenNames,'(Intercept)'));
    pVal2(voxInd) = PVAL(ismember(betweenNames,'(Intercept)'));
end