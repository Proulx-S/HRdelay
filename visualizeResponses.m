function visualizeResponses(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~isfield(p,'chanSpace') || isempty(p.chanSpace)
    p.chanSpace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end

%% Define paths
subjList = p.meta.subjList;
% repoPath = p.paths.repo.in;
%         funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
%             inDir  = 'd';
%             inDir2  = ['e_' p.anaID];
%             inDir3  = 'c';
% %make sure everything is forward slash for mac, linux pc compatibility
% for tmp = {'repoPath' 'funPath' 'inDir'}
%     eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
% end
% clear tmp



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:length(subjList)
    curFile = fullfile(p.dataPath.V1,'resp',[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'resp');
    dAll{subjInd} = resp;
end
d = dAll; clear dAll resp
sessList = fields(d{1});
% Load feature slection
load(fullfile(p.dataPath.V1,'featSel.mat'),'featSel');


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


%% Define feature selection
if ~strcmp(featSel{1,1}.featSeq.info2,p.featSel.global.method)
    error(['p.featSel.global.method and featSel.featSeq.info2 not matching' newline 'Try reruning processFeatSel.m'])
end
featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;


method = 'onlyRetinoFov'; % 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
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
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_act = featSel;


%% Group effect
fTrigGroup = plotTrigGroup(d,p,featSel,1);
fHrGroup = plotHrGroup(d,p,featSel,1);


subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
%% Trigonometric (polar) representation
[fTrig,voxIndTrig,indFovNotAct,indFovAct] = plotTrig(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1);


%% Temporal represeantation
fHr = [];
[fHr,voxIndHr] = plotHr(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1);

%% Time series
curFile = fullfile(p.dataPath.V1,[subjList{p.figOption.subjInd} '.mat']);
fun = load(curFile);
fun = fun.d.fun(p.figOption.sessInd);
runTs = squeeze(cat(5,fun.data{:}));
runLabel = fun.condLabel;
runTs = squeeze(mean(runTs(:,:,runLabel==3),1));
runTs_av = mean(runTs,2)';
runTs_er = bootci(p.boot.n,{@(x)mean(x),runTs'},'Type','percentile');
runTs_er = [runTs_er(2,:) - runTs_av
runTs_av - runTs_er(1,:)];

fTs = figure('windowstyle','docked');
t = 0:119;
hErr = shadedErrorBar(t,runTs_av,runTs_er);
delete(hErr.edge);
ax = gca;
ax.PlotBoxAspectRatio = [120 13 1];
ax.Box = 'off'; ax.Color = 'none';
ax.XTick = 0:12:120;
ax.XGrid = 'on';
ax.XTickLabelRotation = 0;


%% Add-on
figure(fTrig)
ax = findobj(fTrig.Children,'type','Axes');

%info
delay = nan(size(d));
for i = 1:numel(d)
    delay(i) = angle(mean(d{i}.sin(:)));
end
refAngle = delay(p.figOption.subjInd,p.figOption.sessInd);
disp(['ex subj delay=' num2str(-refAngle/pi*6) 's'])
disp(['group mean delay=' num2str(-mean(delay(:))/pi*6) 's'])
disp(['group min delay=' num2str(-min(delay(:))/pi*6) 's'])
disp(['group max delay=' num2str(-max(delay(:))/pi*6) 's'])
[u,v] = pol2cart(refAngle,[0 max(abs(axis))]);
hRef1 = plot(u,v,'k');
[u,v] = pol2cart(wrapToPi([refAngle-pi/2 refAngle+pi/2]),[1 1].*max(abs(axis)));
hRef2 = plot(u,v,'k');

negRatio = nan(size(delay));
for i = 1:numel(d)
    negRatio(i) = nnz( ...
    ( wrapToPi(delay(i)+pi/2) < angle(d{i}.sin(:)) ) & ...
    ( angle(d{i}.sin(:)) < wrapToPi(delay(i)-pi/2) ) ...
    ) ./ numel(d{i}.sin);
end
disp(['ex subj negBOLD%=' num2str(negRatio(subjInd,sessInd).*100) '%'])
disp(['group mean negBOLD%=' num2str(mean(negRatio(:)).*100) '%'])
disp(['group min negBOLD%=' num2str(min(negRatio(:)).*100) '%'])
disp(['group max negBOLD%=' num2str(max(negRatio(:)).*100) '%'])

lg = findobj(fTrig.Children,'type','Legend');
lg.String(end-1:end) = [];

uistack([hRef1 hRef2],'bottom')
hTmp = [findobj(ax.Children,'DisplayName','','type','Patch'); findobj(ax.Children,'DisplayName','','type','Line')];
uistack(hTmp,'bottom')
uistack(findobj(ax.Children,'type','Patch'),'bottom')



figure(fHr)
ax = findobj(fHr.Children,'type','Axes');
hTmp = findobj(ax.Children,'type','Errorbar');
hTmp.Marker = '^';
hTmp.MarkerFaceColor = [1 1 1].*0.6;



%% Save
fullfilename = fullfile(p.figOption.outDir,'Fig5A');
curF = fTrigGroup;
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'eps';
exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); disp([curFile '.' curExt]);
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);

fullfilename = fullfile(p.figOption.outDir,'Fig5B');
curF = fHrGroup;
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'eps';
exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); disp([curFile '.' curExt]);
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);


fullfilename = fullfile(p.figOption.outDir,'Fig2A');
curF = fTrig;
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'eps';
exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); disp([curFile '.' curExt]);
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);

fullfilename = fullfile(p.figOption.outDir,'Fig2B');
curF = fHr;
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'eps';
exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); disp([curFile '.' curExt]);
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);

fullfilename = fullfile(p.figOption.outDir,'Fig1C');
curF = fTs;
curFile = fullfilename;
curExt = 'svg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'eps';
exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); disp([curFile '.' curExt]);
% curF.Color = 'w';
curExt = 'fig';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);
curExt = 'jpg';
saveas(curF,[curFile '.' curExt]); disp([curFile '.' curExt]);


function f = plotTrigGroup(d,p,featSel,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
X = nan([3 size(d)]);
for i = 1:numel(d)
    [x,y,~] = getXYK(d{i},p);
    x = x(:,featSel{i}.indIn);
    for yInd = 1:3
        tmpX = x(y==yInd,:);
        X(yInd,i) = mean(tmpX(:));
    end
end

%% Stats
% tmpX = cat(1,mean(X(1:2,:,:),1),X(3,:,:));
tmpX = angle(X);
tmpX = permute(mean(tmpX,3),[2 1]);
t = table(tmpX(:,1),tmpX(:,2),tmpX(:,3),...
'VariableNames',{'cond1','cond2','cond3'});
Meas = table([1 2 3]','VariableNames',{'cond'});
rm = fitrm(t,'cond1-cond3~1','WithinDesign',Meas);
ranovatbl = ranova(rm);
[H,P,CI,STATS] = ttest(mean(table2array(t(:,1:2)),2),table2array(t(:,3)));
disp(' ')
disp('-------')
disp('Delay')
disp('One-way (3-level) ANOVA')
disp(ranovatbl)
tmp = mean(X(:,:),2);
tmp = (angle(mean(tmp(1:2,:),1)) - angle(tmp(3,:)))/pi*6;
tmp2 = mean(X,3);
tmp2 = (angle(mean(tmp2(1:2,:),1)) - angle(tmp2(3,:)))/pi*6;
disp(['plaid-grat=' num2str(tmp) 'sec (range: ' num2str(min(tmp2),'%0.3f') '-' num2str(max(tmp2),'%0.3f') ')'])
disp(['t=' num2str(STATS.tstat) ', p=' num2str(P)])
disp('-------')
disp(' ')

tmpX = abs(X);
tmpX = permute(mean(tmpX,3),[2 1]);
t = table(tmpX(:,1),tmpX(:,2),tmpX(:,3),...
'VariableNames',{'cond1','cond2','cond3'});
Meas = table([1 2 3]','VariableNames',{'cond'});
rm = fitrm(t,'cond1-cond3~1','WithinDesign',Meas);
ranovatbl = ranova(rm);
[H,P,CI,STATS] = ttest(mean(table2array(t(:,1:2)),2),table2array(t(:,3)));
disp(' ')
disp('-------')
disp('Amplitude')
disp('One-way (3-level) ANOVA on amplitude')
disp(ranovatbl)
tmp = mean(X(:,:),2);
tmp = abs(tmp(3,:)) - abs(mean(tmp(1:2,:),1));
disp(['plaid-grat=' num2str(tmp) '%BOLD'])
disp(['t=' num2str(STATS.tstat) ', p=' num2str(P)])
disp('-------')
disp(' ')

%% Remove random-effect of subject
rho = abs(X) ./ abs(mean(X,1)) .* abs(mean(X(:)));
theta = wrapToPi(angle(X) - angle(mean(X,1)) + angle(mean(X(:))));
[u,v] = pol2cart(theta,rho);
Xnorm = complex(u,v);

%% Average sessions
X = mean(X,3);
Xnorm = mean(Xnorm,3);


if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

Mrkr = 'osd^v><ph';
hPol = cell(1,6);
for subjInd = 1:6
    [u,v] = pol2cart(angle(mean(X(:,subjInd),1)),log(abs(mean(X(:,subjInd),1))+1));
    hPol{yInd,subjInd} = plot(u,v,Mrkr(subjInd),'Color',[1 1 1].*0.5); hold on
    hPol{yInd,subjInd}.MarkerFaceColor = hPol{yInd,subjInd}.Color;
    hPol{yInd,subjInd}.MarkerEdgeColor = 'none';
end

set(gca,'ColorOrderIndex',1)
hPolAv = {};
for yInd = 1:3
    [u,v] = pol2cart(angle(mean(Xnorm(yInd,:),2)),log(abs(mean(Xnorm(yInd,:),2))+1));
    hPolAv{yInd} = plot(u,v,'o'); hold on
    hPolAv{yInd}.MarkerFaceColor = hPolAv{yInd}.Color;
    hPolAv{yInd}.MarkerEdgeColor = 'k';
end

hPolEr = {};    
for yInd = 1:3
    tmpX = permute(Xnorm(yInd,:),[2 1]);
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
    [theta,rho] = cart2pol(polyCont.Vertices(:,1),polyCont.Vertices(:,2));
    [polyCont.Vertices(:,1),polyCont.Vertices(:,2)] = pol2cart(theta,log(rho+1));
    hPolEr{yInd} = plot(polyCont);
    hPolEr{yInd}.LineStyle = 'none';
    hPolEr{yInd}.FaceColor = hPolAv{yInd}.Color;
end

ax = f.Children;

allPoly = findall(ax.Children,'Type','Polygon');
thetaPoly = nan([length(allPoly) 2]);
rhoPoly = nan([length(allPoly) 2]);
xxPoly = nan([length(allPoly) 2]);
yyPoly = nan([length(allPoly) 2]);
for i = 1:length(allPoly)
    xx = allPoly(i).Shape.Vertices(:,1);
    yy = allPoly(i).Shape.Vertices(:,1);
    [theta,rho] = cart2pol(xx,yy);
    rhoPoly(i,:) = [min(rho) max(rho)];
    thetaPoly(i,:) = [min(theta) max(theta)];
    xxPoly(i,:) = [min(xx) max(xx)];
    yyPoly(i,:) = [min(yy) max(yy)];
end
allLine = findall(ax.Children,'Type','Line');
thetaLine = nan([length(allLine) 2]);
rhoLine = nan([length(allLine) 2]);
xxLine = nan([length(allLine) 2]);
yyLine = nan([length(allLine) 2]);
for i = 1:length(allLine)
    xx = allLine(i).XData;
    yy = allLine(i).YData;
    [theta,rho] = cart2pol(xx,yy);
    rhoLine(i,:) = [min(rho) max(rho)];
    thetaLine(i,:) = [min(theta) max(theta)];
    xxLine(i,:) = [min(xx) max(xx)];
    yyLine(i,:) = [min(yy) max(yy)];
end
rhoLim = [min([rhoPoly(:); rhoLine(:)]) max([rhoPoly(:); rhoLine(:)])];
thetaLim = [min([thetaPoly(:); thetaLine(:)]) max([thetaPoly(:); thetaLine(:)])];
xxLim = [min([xxLine(:); xxLine(:)]) max([xxPoly(:); xxLine(:)])];
yyLim = [min([yyLine(:); yyLine(:)]) max([yyPoly(:); yyLine(:)])];


rhoTickMinorVal = [0:0.1:2]';
theta = repmat(linspace(0,2*pi,1000),[length(rhoTickMinorVal) 1]);
[gu,gv] = pol2cart(theta,log(rhoTickMinorVal+1));
hGridRhoMinor = plot(gu',gv','Color',[1 1 1].*0.5);
uistack(hGridRhoMinor,'bottom')
set(hGridRhoMinor,'Color',[1 1 1].*0.8)

rhoTickMajorVal = [0:1:2]';
theta = repmat(linspace(0,2*pi,1000),[length(rhoTickMajorVal) 1]);
[gu,gv] = pol2cart(theta,log(rhoTickMajorVal+1));
hGridRhoMajor = plot(gu',gv','k');
uistack(hGridRhoMajor,'bottom')
uistack(hGridRhoMinor,'bottom')
set(hGridRhoMajor,'Color',[1 1 1].*0.5)

thetaTickMinorVal = 0:0.1:12;
rho = repmat([rhoTickMinorVal(1); rhoTickMinorVal(end)],[1 length(thetaTickMinorVal)]);
% rho = repmat([0; max(axis)],[1 length(thetaTickVal)]);
[gu,gv] = pol2cart(thetaTickMinorVal/12*2*pi,log(rho+1));
hGridThetaMinor = plot(gu,gv,'k');
uistack(hGridThetaMinor,'bottom')
set(hGridThetaMinor,'Color',[1 1 1].*0.8)

thetaTickMajorVal = 0:1:12;
rho = repmat([rhoTickMinorVal(1); rhoTickMinorVal(end)],[1 length(thetaTickMajorVal)]);
% rho = repmat([0; max(axis)],[1 length(thetaTickVal)]);
[gu,gv] = pol2cart(thetaTickMajorVal/12*2*pi,log(rho+1));
hGridThetaMajor = plot(gu,gv,'k');
uistack(hGridThetaMajor,'bottom')
uistack(hGridThetaMinor,'bottom')
set(hGridThetaMajor,'Color',[1 1 1].*0.5)

uistack(allPoly,'top')
uistack(allLine,'top')

xlim(xxLim + [-1 1].*range(xxLim)*0.05)
ylim(yyLim + [-1 1].*range(yyLim)*0.05)

refPts = complex(mean(xxLim),mean(yyLim));
[theta,rho] = meshgrid(thetaTickMinorVal,rhoTickMinorVal);
[xx,yy] = pol2cart(theta/12*2*pi,log(rho+1));
TickMinorVal = complex(xx(:),yy(:));
[a,b] = sort(abs(TickMinorVal - refPts),'ascend');
for i = 1:4
%     plot(TickMinorVal(b(i)),'or');
    tmpText = [num2str(abs(angle(TickMinorVal(b(i)))/pi*6)) 'sec, ' num2str(exp(abs(TickMinorVal(b(i))))-1) '%BOLD'];
    text(real(TickMinorVal(b(i))),imag(TickMinorVal(b(i))),tmpText)
end

ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

hPol = [hPol{:}]';
hPolAv = [hPolAv{:}]';
hPolEr = [hPolEr{:}]';
hGridThetaMajor;
hGridThetaMinor;
hGridRhoMajor;
hGridRhoMinor;

uistack(hPolEr)
uistack(hPol)
uistack(hPolAv)
set(hPolAv,'MarkerSize',5)
set(hPolAv,'MarkerEdgeColor','none')


function f = plotHrGroup(d,p,featSel,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
X = nan([3 size(d)]);
for i = 1:numel(d)
    [x,y,~] = getXYK(d{i},p);
    x = x(:,featSel{i}.indIn);
    for yInd = 1:3
        tmpX = x(y==yInd,:);
        X(yInd,i) = mean(tmpX(:));
    end
end
x = X; clear X
x = permute(x,[1 4 2 3]);

%% Remove random-effect of subject
rho = abs(x) ./ abs(mean(x,1)) .* abs(mean(x(:)));
theta = wrapToPi(angle(x) - angle(mean(x,1)) + angle(mean(x(:))));
[u,v] = pol2cart(-theta,rho);
x = complex(u,v);

%% Plot fits
tmpX = squeeze(mean(x,4))';
t = 0:0.1:12;
tmpXsin = nan([length(t) size(tmpX)]);
tmpXsin(:,:) = ( real(tmpX(:))*sin(t/6*pi) ...
    + imag(tmpX(:))*cos(t/6*pi) )';
% tmpXsin(:,:) = ( real(tmpX(:))*sin(linspace(0,2*pi,13)) ...
%     + imag(tmpX(:))*cos(linspace(0,2*pi,13)) )';
tmpXsin = permute(tmpXsin,[2 3 1]);

%     tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
%     tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
n = size(tmpX,1);
tmpXsinMean = mean(tmpXsin(:,:),1);
tmpXsinBoot = nan(nBoot,size(tmpXsin(:,:),2));
for bootInd = 1:nBoot
    boot = nan(n,1);
    for i = 1:n
        boot(i) = randperm(n,1);
    end
    tmpXsinBoot(bootInd,:) = mean(tmpXsin(boot,:),1);
end
tmpXsinEr = prctile(tmpXsinBoot,[97.5 2.5],1);
tmpXsinEr(1,:) = tmpXsinEr(1,:) - tmpXsinMean;
tmpXsinEr(2,:) = tmpXsinMean - tmpXsinEr(2,:);

hrFitAv = nan([3 length(t)]);
hrFitEr = nan([2 3 length(t)]);
hrFitAv(:) = tmpXsinMean;
hrFitEr(:,:) = tmpXsinEr;

if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end
hTmp = plot([1 2 3; 1 2 3]);
cList = cat(1,hTmp.Color);
delete(hTmp)
hShaded = {};
for yInd = 1:3
    hShaded{yInd} = shadedErrorBar(t,hrFitAv(yInd,:),squeeze(hrFitEr(:,yInd,:)),'lineprops',{'Color' cList(yInd,:)}); hold on
    delete(hShaded{yInd}.edge)
end
%     hPlot2vox{yInd} = plot(t,tmpXmean,'Color',cList(yInd,:)); hold on
%     hPlot1vox{yInd}
%     hPlot1vox{yInd} = plot(t,tmpX,'-','Color',cList(yInd,:)); hold on
hShaded = [hShaded{:}];

%% Plot individual subjects
hShadedSubj = {};
for i = 1:size(d,1)
    hr1 = permute(squeeze(mean(d{i,1}.hr(featSel{i}.indIn,:,:,:,:,:),1)),[3 1 2]);
    hr1 = hr1 - mean(hr1,1);
    hr2 = permute(squeeze(mean(d{i,2}.hr(featSel{i}.indIn,:,:,:,:,:),1)),[3 1 2]);
    hr2 = hr2 - mean(hr2,1);
    hr = cat(2,hr1,hr2);
    hr = hr(:,:);
%     hrAv(:,i) = mean(hr,2);
    hrAv = mean(hr,2);
    hrEr = nan(size(hrAv,1),2);
    for ii = 1:size(hr,1)
        hrEr(ii,:) = bootci(p.boot.n,{@mean,hr(ii,:)},'Type','per','Alpha',0.05);
    end
    hrEr(:,1) = hrEr(:,1) - hrAv;
    hrEr(:,2) = hrAv - hrEr(:,2);
    
    hrAv = cat(1,hrAv,hrAv(1));
    hrEr = cat(1,hrEr,hrEr(1,:));
    
    t2 = (1:length(hrAv))-1;
    
    hShadedSubj{i} = shadedErrorBar(t2',hrAv,hrEr,'lineprops',{'Color' [1 1 1].*0.8}); hold on
    delete(hShadedSubj{i}.edge)
end
hShadedSubj = [hShadedSubj{:}];
mainLine = [hShadedSubj.mainLine];
patch = [hShadedSubj.patch];
uistack(mainLine,'bottom')
uistack(patch,'bottom')

%% Decorations
set(mainLine,'Color',[1 1 1].*0.7)
set(patch,'FaceColor',[1 1 1].*0.7)
set(patch,'Visible','on')
ax = gca;
ax.PlotBoxAspectRatio = [1.1 1 1];
ax.XTick = 0:3:12;
ax.YTick = -3:1:3;
ax.YLim(2) = -ax.YLim(1);
grid on
ax.Box = 'off';



function [f,voxInd,indFovNotAct,indFovAct] = plotTrig(d,p,featSel1,featSel2,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
[X,y,~] = getXYK(d,p);

% [u,v] = pol2cart(angle(X),log(abs(X)+1));
[u,v] = pol2cart(angle(X),abs(X));
x = complex(u,v);

indFovNotAct = featSel2.indIn & ~featSel1.indIn;
indFovAct = featSel1.indIn;
xFovNotAct = x(:,indFovNotAct);
xFovAct = x(:,indFovAct);


if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

[u,v] = pol2cart(angle(mean(xFovNotAct,1)),log(abs(mean(xFovNotAct,1))+1));
hPolFov = plot(u,v,'.'); hold on
% hPolFov = polar(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
% hPolFov = polarplot(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
hPolFov.MarkerFaceColor = 'k';
hPolFov.MarkerEdgeColor = 'k';
hPolFov.MarkerSize = eps;


[u,v] = pol2cart(angle(mean(xFovAct,1)),log(abs(mean(xFovAct,1))+1));
hPolAct= plot(u,v,'.'); hold on
% hPolAct = polar(angle(mean(xAct,1)),abs(mean(xAct,1)),'.'); hold on
% hPolAct = polarplot(angle(mean(xAct,1)),abs(mean(xAct,1)),'.'); hold on
hPolAct.MarkerFaceColor = [1 1 1]-eps;
hPolAct.MarkerEdgeColor = [1 1 1]-eps;
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
voxInd = false(size(act));
voxInd(b) = true;


ax = gca;
ax.ColorOrderIndex = 1;

yList = unique(y);
hPol1vox = {};
for yInd = 1:length(yList)
    tmpX = squeeze(x(yList(yInd)==y,b));
    tmpXmean = mean(tmpX,1);
    [u,v] = pol2cart(angle(tmpXmean),log(abs(tmpXmean)+1));
    hPol1vox{yInd} = plot(u,v,'.'); hold on
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
    [theta,rho] = cart2pol(polyCont.Vertices(:,1),polyCont.Vertices(:,2));
    [polyCont.Vertices(:,1),polyCont.Vertices(:,2)] = pol2cart(theta,log(rho+1));
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

[~,rho] = cart2pol(hPolFov.XData,hPolFov.YData);
ind = rho>log(rhoTickVal(end)+1);
hPolFov.XData(ind) = []; hPolFov.YData(ind) = [];

theta = repmat(linspace(0,2*pi,100),[length(rhoTickVal) 1]);
[gu,gv] = pol2cart(theta,log(rhoTickVal+1));
hGridConc = plot(gu',gv','k');


thetaTickVal = 0:1:11;
rho = repmat([rhoTickVal(1); rhoTickVal(end)],[1 length(thetaTickVal)]);
% rho = repmat([0; max(axis)],[1 length(thetaTickVal)]);
[gu,gv] = pol2cart(thetaTickVal/12*2*pi,log(rho+1));
hGridRad1 = plot(gu,gv,'k');

thetaTickVal = 0:3:9;
rho = repmat([0; rhoTickVal(1)],[1 length(thetaTickVal)]);
% rho = repmat([0; max(axis)],[1 length(thetaTickVal)]);
[gu,gv] = pol2cart(thetaTickVal/12*2*pi,log(rho+1));
hGridRad2 = plot(gu,gv,'k');
hGridRad = [hGridRad1; hGridRad2];


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




function [f,voxInd] = plotHr(d,p,featSel1,featSel2,visibilityFlag)
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
voxInd = false(size(act));
voxInd(b) = true;


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


%% Decorations
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



