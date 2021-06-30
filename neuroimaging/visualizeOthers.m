function visualizeOthers(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
% if ~exist('figOption','var') || isempty(figOption)
%     figOption.save = 0;
%     figOption.subj = 1; % 'all' or subjInd
% end
if ~isfield(p,'chanSpace') || isempty(p.chanSpace)
    p.chanSpace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end
% tic

%% Define paths
subjList = p.meta.subjList;
repoPath = p.paths.repo.in;
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            inDir2  = ['e_' p.anaID];
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:length(subjList)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
    dAll{subjInd} = res;
end
d = dAll; clear dAll
sessList = fields(d{1});
% Load feature slection
load(fullfile(funPath,inDir2,'featSel.mat'),'featSel');
% if verbose
%     disp('---');
%     disp(['Channel space: ' p.chanSpace '-' p.condPair]);
%     disp(['Complex space: ' p.complexSpace]);
% end

%% Reorganize
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
        dP{subjInd,sessInd}.featSel = featSel{subjInd,sessInd};
    end
end
d = dP; clear dP

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Define feature selection
if ~strcmp(featSel{1,1}.featSeq.info2,p.featSel.global.method)
    error(['p.featSel.global.method and featSel.featSeq.info2 not matching' newline 'Try reruning processFeatSel.m'])
end
featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;


method = 'onlyRetinoFov'; % 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
% [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,p.featSel.global.method,p.condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_fov = featSel;


method = 'upToActivation'; % 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
% [ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,p.featSel.global.method,p.condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_act = featSel;

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc



subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
%% Trigonometric (polar) representation
fTrig = [];
if p.figOption.verbose==1
    fTrig = [fTrig plotTrig(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},0)];
elseif p.figOption.verbose>1
    fTrig = [fTrig plotTrig(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1)];
end

%% Temporal represeantation
fHr = [];
if p.figOption.verbose==1
    fHr = [fHr plotHr(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},0)];
elseif p.figOption.verbose>1
    fHr = [fHr plotHr(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1)];
end

%% Add-on
% figure(fTrig)
% ax = findobj(fTrig.Children,'type','PolarAxes');
% tmp = ax.Children([4 5]);
% theta = [tmp(1).ThetaData tmp(2).ThetaData];
% rho = [tmp(1).RData tmp(2).RData];
% [u,v] = pol2cart(theta,rho);
% tmp = complex(u,v);
% tmp = mean(tmp);
% hRef(1) = polarplot([0 angle(tmp)],ax.RLim,'k-');
% theta = wrapToPi([-1 1].*pi/2 + angle(tmp));
% rho = [1 1].*ax.RLim(2);
% hRef(2) = polarplot(theta,rho,'k-');




%% Save
fullfilename = fullfile(p.figOption.finalDir,'polarRespVec');
curF = fTrig;
% curA = findobj(curF.Children,'type','PolarAxes');
% axColor = curA.Color;
% curF.Color = 'none';
% curA.Color = axColor;
% set(findobj(curF.Children,'type','Axes'),'color','none')
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
% curF.Color = 'w';
% curA.Color = axColor;
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end

fullfilename = fullfile(p.figOption.finalDir,'Hr');
curF = fHr;
% curF.Color = 'none';
% set(findobj(curF.Children,'type','Axes'),'color','none')
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
% curF.Color = 'w';
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end


function f = plotTrig(d,p,featSel1,featSel2,visibilityFlag)
nBoot = p.boot.n;
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

hPolFov = plot(real(mean(xFov,1)),imag(mean(xFov,1)),'.'); hold on
% hPolFov = polar(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
% hPolFov = polarplot(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
hPolFov.MarkerFaceColor = 'k';
hPolFov.MarkerEdgeColor = 'k';
% hPolFov.MarkerFaceColor = [1 1 1].*0.5;
% hPolFov.MarkerEdgeColor = [1 1 1].*0.5;
hPolFov.MarkerSize = eps;

hPolAct = plot(real(mean(xAct,1)),imag(mean(xAct,1)),'.'); hold on
% hPolAct = polar(angle(mean(xAct,1)),abs(mean(xAct,1)),'.'); hold on
% hPolAct = polarplot(angle(mean(xAct,1)),abs(mean(xAct,1)),'.'); hold on
hPolAct.MarkerFaceColor = [1 1 1]-eps;
hPolAct.MarkerEdgeColor = [1 1 1]-eps;
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
    tmpX = squeeze(x(yList(yInd)==y,b));
    tmpXmean = mean(tmpX,1);
    hPol1vox{yInd} = plot(real(tmpXmean),imag(tmpXmean),'.'); hold on
%     hPol1vox{yInd} = polar(angle(tmpXmean),abs(tmpXmean),'.'); hold on
    %     hPol1vox{yInd} = polarplot(angle(mean(x(yList(yInd)==y,b),1)),abs(mean(x(yList(yInd)==y,b),1)),'.'); hold on
end

hPol2vox = {};    
for yInd = 1:length(yList)
    tmpX = squeeze(x(yList(yInd)==y,b));
    n = size(tmpX,1);
    tmpXboot = nan(nBoot,size(tmpX,2));
    for bootInd = 1:nBoot
        boot = nan(n,1);
        for i = 1:n
            boot(i) = randperm(n,1);
        end
        tmpXboot(bootInd,:) = mean(tmpX(boot,:),1);
    end
    polyCont = credibleInt2D([real(tmpXboot) imag(tmpXboot)],0.05);
%     hPol1vox{yInd} = polar(angle(mean(x(yList(yInd)==y,b),1)),abs(mean(x(yList(yInd)==y,b),1)),'.'); hold on
    hPol2vox{yInd} = plot(polyCont);
    hPol2vox{yInd}.LineStyle = 'none';
    hPol2vox{yInd}.FaceColor = hPol1vox{yInd}.Color;    
end
set([hPol1vox{:}],'MarkerSize',15)


ax.PlotBoxAspectRatio = [1 1 1];
ax.DataAspectRatio = [1 1 1];

xLim = xlim;
yLim = ylim;


% rhoTickVal = [0:0.1:0.9 1:9 10:10:100]';
rhoTickVal = [1:9 10]';
theta = repmat(linspace(0,2*pi,100),[length(rhoTickVal) 1]);
[gu,gv] = pol2cart(theta,log(rhoTickVal+1));
hGridConc = plot(gu',gv','k');

thetaTickVal = 0:1:11;
rho = repmat([0; rhoTickVal(end)],[1 length(thetaTickVal)]);
% rho = repmat([0; max(axis)],[1 length(thetaTickVal)]);
[gu,gv] = pol2cart(thetaTickVal/12*2*pi,log(rho+1));
hGridRad = plot(gu,gv,'k');

[gu,gv] = pol2cart(theta(1,:),log(rhoTickVal(end)+1));
hBck = patch(gu,gv,[1 1 1].*0.6);
uistack(hGridRad,'bottom')
uistack(hGridConc,'bottom')
uistack(hBck,'bottom')
uistack(flip([hPol2vox{:}],2),'top')
uistack([hPol1vox{:}],'top')
hPol1vox = [hPol1vox{:}];
for i = 1:length(hPol2vox)
    hPol1vox(i).Marker = 'o';
    hPol1vox(i).MarkerFaceColor = hPol2vox{i}.FaceColor;
    hPol1vox(i).MarkerEdgeColor = 'k';
    hPol1vox(i).MarkerSize = 5;
    hPol1vox(i).LineWidth = eps;
end

set(hGridRad,'LineWidth',eps);
set(hGridRad,'Color',[1 1 1].*0.5);
set(hGridConc,'LineWidth',eps);
set(hGridConc,'Color',[1 1 1].*0.5);

xlim([-1 1].*log(rhoTickVal(end)+1))
ylim([-1 1].*log(rhoTickVal(end)+1))
ax.XTick = [];
ax.YTick = [];
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

ax.Color = 'none';
f.Color = 'w';

legend([hPolFov hPolAct hPol1vox],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid'},'box','off')
% box off
% % ylim([-1 1].*max(abs([xLim yLim])))
% % xlim([-1 1].*max(abs([xLim yLim])))
% 
% rhoTickVal = [0:0.1:0.9 1:9 10:10:100];
% ax
% ax.RAxis.TickValues = log(rhoTickVal+1);
% 
% tickPos = ax.RAxis.TickValues(1:length(ax.RAxis.TickLabels));
% rhoTickVal = round(exp(tickPos)-1,1,'significant');
% 
% labelVal = [0 1 10];
% ind = ismember(rhoTickVal,labelVal);
% labelVal = cellstr(num2str(rhoTickVal'))';
% labelVal(~ind) = {''};
% ax.RAxis.TickLabels = labelVal;
% 
% 
% theta = 0:0.01:2*pi;
% rho = ones(size(theta)).*log(1+1);
% hPol1 = polarplot(theta,rho,'k');
% rho = ones(size(theta)).*log(10+1);
% hPol10 = polarplot(theta,rho,'k');
% uistack([hPol1 hPol10],'bottom')
% 
% legend([hPolFov hPolAct hPol1vox{:}],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid'},'box','off')
% % f.Color = 'w';
% ax.Box = 'off';
% 
% ax.Color = [1 1 1].*0.6;






function f = plotHr(d,p,featSel1,featSel2,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
[xSin,~,~] = getXYK(d,p,0);
[x,y,~] = getXYK(d,p,1);
x = cat(3,x,x(:,:,1));
xFov = x(:,featSel2.indIn & ~featSel1.indIn,:);
xAct = x(:,featSel1.indIn,:);
t = (1:size(xAct,3))-1;

if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

tmpX = squeeze(mean(xFov,1));
% tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
% tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
hPlotFov = plot(t,tmpX,'Color','k'); hold on
set(hPlotFov,'Color',[get(hPlotFov(1),'Color') 0.15])
set(hPlotFov,'LineWidth',eps);

tmpX = squeeze(mean(xAct,1));
% tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
% tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
hPlotAct = plot(t,tmpX,'w'); hold on
set(hPlotAct,'Color',[get(hPlotAct(1),'Color') 0.15])
hPolAct.LineWidth = eps;

ax = gca;
ax.Color = [1 1 1].*0.6;
axis tight


%% Choose one voxel
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



%% Plot the one voxel
tmpX = squeeze(x(:,b,:));
n = size(tmpX,1);
tmpXmean = mean(tmpX,1);
tmpXboot = nan(nBoot,size(tmpX,2));
for bootInd = 1:nBoot
    boot = nan(n,1);
    for i = 1:n
        boot(i) = randperm(n,1);
    end
    tmpXboot(bootInd,:) = mean(tmpX(boot,:),1);
end
tmpXer = prctile(tmpXboot,[97.5 2.5],1);
tmpXer(1,:) = tmpXer(1,:) - tmpXmean;
tmpXer(2,:) = tmpXmean - tmpXer(2,:);
hPlot2vox = errorbar(t,tmpXmean,tmpXmean-tmpXer(2,:),tmpXer(1,:)-tmpXmean,'ok');
hPlot2vox.CapSize = 0;
hPlot2vox.MarkerEdgeColor = 'k';
hPlot2vox.MarkerFaceColor = 'w';%hPlot2vox.MarkerEdgeColor;



%% Plot conditions for the one voxel
yList = unique(y);
hPlot1vox = {};
cList = colororder;
for yInd = 1:length(yList)
    tmpX = squeeze(x(yList(yInd)==y,b,:));
    tmpXsin = squeeze(xSin(yList(yInd)==y,b));
    tmpXsin = real(tmpXsin)*sin(linspace(0,2*pi,13)) ...
        + imag(tmpXsin)*cos(linspace(0,2*pi,13));
    
%     tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
%     tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
    n = size(tmpX,1);
    tmpXmean = mean(tmpX,1);
    tmpXsinMean = mean(tmpXsin,1);
    tmpXboot = nan(nBoot,size(tmpX,2));
    tmpXsinBoot = nan(nBoot,size(tmpXsin,2));
    for bootInd = 1:nBoot
        boot = nan(n,1);
        for i = 1:n
            boot(i) = randperm(n,1);
        end
        tmpXboot(bootInd,:) = mean(tmpX(boot,:),1);
        tmpXsinBoot(bootInd,:) = mean(tmpXsin(boot,:),1);
    end
    tmpXer = prctile(tmpXboot,[97.5 2.5],1);
    tmpXer(1,:) = tmpXer(1,:) - tmpXmean;
    tmpXer(2,:) = tmpXmean - tmpXer(2,:);
    tmpXsinEr = prctile(tmpXsinBoot,[97.5 2.5],1);
    tmpXsinEr(1,:) = tmpXsinEr(1,:) - tmpXsinMean;
    tmpXsinEr(2,:) = tmpXsinMean - tmpXsinEr(2,:);
    hPlot1vox{yInd} = shadedErrorBar(t,tmpXsinMean,tmpXsinEr,'lineprops',{'Color' cList(yInd,:)}); hold on
    delete(hPlot1vox{yInd}.edge)
%     hPlot2vox{yInd} = plot(t,tmpXmean,'Color',cList(yInd,:)); hold on
%     hPlot1vox{yInd}
%     hPlot1vox{yInd} = plot(t,tmpX,'-','Color',cList(yInd,:)); hold on
end
hPlot1vox = [hPlot1vox{:}];


%% Decortations
tmpX = mean(xAct,1);
ylim([-1 1].*max(abs(tmpX(:))))
% ylim([-1/3 1/3].*max(abs(ylim)))

% uistack(hPlot1vox,'top')
legend([hPlotFov(1) hPlotAct(1) [hPlot1vox.mainLine] [hPlot1vox.patch]],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid' '95' '95' '95'},'box','off','Location','NorthWest')
% legend([hPlotFov(1) hPlotAct(1) hPlot1vox(1,:)],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid'},'box','off')
% f.Color = 'w';
ax = gca;
ax.Box = 'off';

% ax.XAxis.TickValues = 0:1:11;
ax.XAxis.TickValues = 0:3:12;



