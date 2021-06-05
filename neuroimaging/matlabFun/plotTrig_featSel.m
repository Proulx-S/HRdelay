function f = plotTrig_featSel(d,p,featSel1,featSel2,visibilityFlag)
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
[X,y,~] = getXYK(d,p);

[u,v] = pol2cart(angle(X),log(abs(X)+1));
x = complex(u,v);

xFov = x(:,featSel2.indIn & ~featSel1.indIn);
xAct = x(:,featSel1.indIn);

if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

hPolFov = polarplot(angle(mean(xFov,1)),abs(mean(xFov,1)),'.k'); hold on
hPolFov.MarkerFaceColor = [1 1 1].*0.5;
hPolFov.MarkerEdgeColor = [1 1 1].*0.5;
hPolFov.MarkerSize = eps;

hPolAct = polarplot(angle(mean(xAct,1)),abs(mean(xAct,1)),'.k'); hold on
% hPolAct.MarkerFaceColor = [1 1 1].*0.5
% hPolAct.MarkerEdgeColor = [1 1 1].*0.5
hPolAct.MarkerSize = eps;


featSelSteps_labelList = featSel1.featSeq.featSelList;
for i = 1:length(featSelSteps_labelList)
    tmp = strsplit(featSelSteps_labelList{i},':');
    featSelSteps_labelList(i) = tmp(1);
end

featSelStep_ind = ismember(featSelSteps_labelList,{'act' 'respVecSig'});
condPair_ind = logical([1 0 0 0]);
if ~all(featSel1.featSeq.condPairList{condPair_ind}==[1 2 3])
    error('lazy coding')
end

act = sqrt(sum(featSel1.featSeq.featVal(:,featSelStep_ind,condPair_ind).^2,2));
[~,b] = min(abs(prctile(act(featSel1.indIn),80)-act));


ax = gca;
ax.ColorOrderIndex = 1;

yList = unique(y);
hPol1vox = {};
for yInd = 1:length(yList)
    hPol1vox{yInd} = polarplot(angle(x(yList(yInd)==y,b)),abs(x(yList(yInd)==y,b)),'.'); hold on
end
set([hPol1vox{:}],'MarkerSize',15)

tickVal = [0:0.1:0.9 1:9 10:10:100];
ax.RAxis.TickValues = log(tickVal+1);

tickPos = ax.RAxis.TickValues(1:length(ax.RAxis.TickLabels));
tickVal = round(exp(tickPos)-1,1,'significant');

labelVal = [0 1 10];
ind = ismember(tickVal,labelVal);
labelVal = cellstr(num2str(tickVal'))';
labelVal(~ind) = {''};
ax.RAxis.TickLabels = labelVal;


theta = 0:0.01:2*pi;
rho = ones(size(theta)).*log(1+1);
hPol1 = polarplot(theta,rho,'k');
rho = ones(size(theta)).*log(10+1);
hPol10 = polarplot(theta,rho,'k');
uistack([hPol1 hPol10],'bottom')

legend([hPolFov hPolAct hPol1vox{:}],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid'},'box','off')
f.Color = 'w';
ax.Box = 'off';



