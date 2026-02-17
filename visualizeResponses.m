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


% No feature selection: Include all V1 ROI voxels
method = 'v1ROI'; % 'v1ROI' 'onlyRetinoFov' 'upToActivation'
condPair = 'all';
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,method,condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end
featSel_v1 = featSel;
featSel_v1{1,1}.featSeq.featSelList(ind_nSpecFeatSel)

% Minimal feature selection: Include all V1 ROI voxels then select voxels (1) within stimulus FOV
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
featSel_fov{1,1}.featSeq.featSelList(ind_nSpecFeatSel)


% Final feature selection: Include all V1 ROI voxels then select voxels (1) within stimulus FOV, (2) that are significantly activated, (3) least likely to be veins and (4) most discriminative
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
featSel_act{1,1}.featSeq.featSelList(ind_nSpecFeatSel)


%% Group effect
[fTrigGroup,fSuppCorrelation] = plotTrigGroup(d,p,featSel,1); % Fig5A, FigSuppCorrelation
fHrGroup = plotHrGroup(d,p,featSel,1); % Fig5B


subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
%% Trigonometric (polar) representation
% [fTrig,voxIndTrig,indFovNotAct,indFovAct] = plotTrig(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1);
[fTrig,voxIndTrig,indFovNotAct,indFovAct] = plotTrig(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_v1{subjInd,sessInd},1);


%% Temporal represeantation
% fHr = [];
% [fHr,voxIndHr]       = plotHr(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_fov{subjInd,sessInd},1);
fHr = [];
[fHr,voxIndHr] = plotHr(d{subjInd,sessInd},p,featSel_act{subjInd,sessInd},featSel_v1{subjInd,sessInd},voxIndTrig,1);

%% Time series
curFile = fullfile(p.dataPath.V1,'ts',[subjList{p.figOption.subjInd} '.mat']);
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
% hTmp.Marker = '^';
% hTmp.MarkerFaceColor = [1 1 1].*0.6;



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

fullfilename = fullfile(p.figOption.outDir,'FigSuppCorrelation');
curF = fSuppCorrelation;
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


function [f,f2] = plotTrigGroup(d,p,featSel,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
X    = nan([3 size(d)]);
Xrun = cell(size(d));
for i = 1:numel(d)
    [x,y,~] = getXYK(d{i},p);
    x = x(:,featSel{i}.indIn);
    % sort conditions
    for yInd = 1:3
        tmpX = x(y==yInd,:);
        X(yInd,i) = mean(tmpX(:));
        Xrun{i}(yInd,:) = mean(tmpX,2);
    end
end
for subjInd = 1:size(d,1)
    Xrun{subjInd,1} = cat(2,Xrun{subjInd,:});
end
Xrun(:,2) = [];

%% Delay estimate precision
phaseWrap = nan(size(d,1),3);
delayCI = nan(size(d,1),3,2);
for subjInd = 1:size(d,1)
    for condInd = 1:3
        phaseWrap(subjInd,condInd) = any(angle(Xrun{subjInd}(condInd,:))>0);
        delayCI(subjInd,condInd,:) = bootci(p.boot.n,{@mean,-angle(Xrun{subjInd}(condInd,:))./pi*6},'Type','per');
    end
end
any(phaseWrap(:)); % no phase wrap
delayCIwidth = diff(delayCI,[],3);
disp(['delay CI width (mean)=' num2str(mean(delayCIwidth(:))) 'sec'])
disp(['delay CI width (min)=' num2str(min(delayCIwidth(:))) 'sec'])
disp(['delay CI width (max)=' num2str(max(delayCIwidth(:))) 'sec'])

%% Stats -- main effects
disp(' ')
disp('-------')
disp('Delay')
% tmpX = angle(X); % take the delay
% tmpX = permute(mean(tmpX,3),[2 1]); % average across sessions (no phase wrap)
tmpX = permute(mean(X,3),[2 1]); % average across sessions
tmpX = angle(tmpX); % take the delay
disp('One-way (3-level) ANOVA')
t = table(tmpX(:,1),tmpX(:,2),tmpX(:,3),...
'VariableNames',{'cond1','cond2','cond3'});
Meas = table([1 2 3]','VariableNames',{'cond'});
rm = fitrm(t,'cond1-cond3~1','WithinDesign',Meas);
ranovatbl = ranova(rm);
disp(ranovatbl)
delayGrat   = -mean(table2array(t(:,1:2)),2)./pi*6; % average orientations and convert to seconds (no phase wrap)
delayPlaid = -table2array(t(:,3))./pi*6; % convert to seconds (no phase wrap)
[H,pTtest,CI,STATS] = ttest(delayPlaid,delayGrat);
disp(['plaid-grat=' num2str(mean(delayPlaid - delayGrat)) 'sec (range: ' num2str(min(delayPlaid - delayGrat),'%0.3f') ' to ' num2str(max(delayPlaid - delayGrat),'%0.3f') ')'])
disp(['t=' num2str(STATS.tstat) ', p=' num2str(pTtest)])
[pSignedRank, h, stats] = signrank(delayPlaid,delayGrat);
disp(['signedRank=' num2str(stats.signedrank) ', p=' num2str(pSignedRank)])
disp('-------')
disp(' ')

disp(' ')
disp('-------')
disp('Amplitude')
% tmpX = abs(X); % take the amplitude
% tmpX = permute(mean(tmpX,3),[2 1]); % average across sessions
tmpX = permute(mean(X,3),[2 1]); % average across sessions
tmpX = abs(tmpX); % take the amplitude
disp('One-way (3-level) ANOVA')
t = table(tmpX(:,1),tmpX(:,2),tmpX(:,3),...
'VariableNames',{'cond1','cond2','cond3'});
Meas = table([1 2 3]','VariableNames',{'cond'});
rm = fitrm(t,'cond1-cond3~1','WithinDesign',Meas);
ranovatbl = ranova(rm);
disp(ranovatbl)
ampGrat  = mean(table2array(t(:,1:2)),2); % average orientations and convert to seconds (no phase wrap)
ampPlaid = table2array(t(:,3)); % convert to seconds (no phase wrap)
[H,pTtest,CI,STATS] = ttest(ampPlaid,ampGrat);
disp(['plaid-grat=' num2str(mean(ampPlaid - ampGrat)) '%BOLD (range: ' num2str(min(ampPlaid - ampGrat),'%0.3f') ' to ' num2str(max(ampPlaid - ampGrat),'%0.3f') ')'])
disp(['t=' num2str(STATS.tstat) ', p=' num2str(pTtest)])
[pSignedRank, h, stats] = signrank(ampPlaid,ampGrat);
disp(['signedRank=' num2str(stats.signedrank) ', p=' num2str(pSignedRank)])
disp('-------')
disp(' ')


%% Plot -- main effects

% Remove random-effect of subject
rho = abs(X) ./ abs(mean(X,1)) .* abs(mean(X(:)));
theta = wrapToPi(angle(X) - angle(mean(X,1)) + angle(mean(X(:))));
[u,v] = pol2cart(theta,rho);
Xnorm = complex(u,v);

% Average sessions
X     = mean(X,3);
Xnorm = mean(Xnorm,3);


% Plot
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
[gu,gv] = pol2cart(thetaTickMinorVal/12*2*pi,log(rho+1));
hGridThetaMinor = plot(gu,gv,'k');
uistack(hGridThetaMinor,'bottom')
set(hGridThetaMinor,'Color',[1 1 1].*0.8)

thetaTickMajorVal = 0:1:12;
rho = repmat([rhoTickMinorVal(1); rhoTickMinorVal(end)],[1 length(thetaTickMajorVal)]);
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


%% Stats and plot -- delay vs amp correlation (Supplementary Figure 3)
for subjInd = 1:6
    delayPlaidRun{subjInd} = -                    angle(Xrun{subjInd}(3  ,:))       ./pi*6; % no phase wrap
    delayGratRun{subjInd}  = -angle( mean( exp(1i*angle(Xrun{subjInd}(1:2,:))) ,1) )./pi*6; % no phase wrap
    ampPlaidRun{subjInd}   =                        abs(Xrun{subjInd}(3  ,:))       ;
    ampGratRun{subjInd}    =         mean(          abs(Xrun{subjInd}(1:2,:))  ,1  );    
end
for subjInd = 1:6
    delayEffect{subjInd} = delayPlaidRun{subjInd} - delayGratRun{subjInd};
    ampEffect{subjInd}   = ampPlaidRun{subjInd}   - ampGratRun{subjInd};
end
for subjInd = 1:6
    delayEffectAv(subjInd) = mean(delayEffect{subjInd});
    er = bootci(p.boot.n,{@mean,delayEffect{subjInd}},'Type','percentile');
    delayEffectErNeg(subjInd) = delayEffectAv(subjInd) - er(1,:);
    delayEffectErPos(subjInd) = er(2,:) - delayEffectAv(subjInd);
    ampEffectAv(subjInd) = mean(ampEffect{subjInd});
    er = bootci(p.boot.n,{@mean,ampEffect{subjInd}},'Type','percentile');
    ampEffectErNeg(subjInd) = ampEffectAv(subjInd) - er(1,:);
    ampEffectErPos(subjInd) = er(2,:) - ampEffectAv(subjInd);
end

disp(' ')
disp('-------')
disp('Delay vs Amplitude (averaging runs instead of sessions)')
[rPearson,pPearson]   = corr(delayEffectAv',ampEffectAv','Type','Pearson');
[rSpearman,pSpearman] = corr(delayEffectAv',ampEffectAv','Type','Spearman');
disp(['Pearson''s rho=' num2str(rPearson,'%0.2f') ', p=' num2str(pPearson,'%0.2f')])
disp(['Spearman''s rho=' num2str(rSpearman,'%0.2f') ', p=' num2str(pSpearman,'%0.2f')])
disp('-------')
disp(' ')


f2 = figure;
hScat = scatter(delayEffectAv,ampEffectAv); hold on
hScat.Visible = 'off';
hErr = errorbar(delayEffectAv,ampEffectAv,delayEffectErNeg,delayEffectErPos,ampEffectErNeg,ampEffectErPos);
hErr.CapSize = 0;
hErr.LineStyle = 'none';
hErr.Marker = 'o';
hErr.MarkerEdgeColor = 'none';
hErr.Color = 'k';
hErr.MarkerFaceColor = 'k';
grid on
grid minor
axis tight square
hLine = lsline;
hLine.LineStyle = '--';
xLim = xlim; xLim = [-1.05 1.05].*max(abs(xLim)); xlim(xLim)
yLim = ylim; yLim = [-1.05 1.05].*max(abs(yLim)); ylim(yLim)
xline(0,'k'); yline(0,'k');
xlabel('plaid-grat delay difference (s)')
ylabel('plaid-grat amplitude difference (%BOLD)')
title({
    'Amplitude vs. Delay (averaging runs instead of sessions)'
    ['Pearson''s rho=' num2str(rPearson,'%0.2f') ', p=' num2str(pPearson,'%0.2f')]
    ['Spearman''s rho=' num2str(rSpearman,'%0.2f') ', p=' num2str(pSpearman,'%0.2f')]
})







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
hPolFov = plot(u,v,'o'); hold on
% hPolFov = polar(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
% hPolFov = polarplot(angle(mean(xFov,1)),abs(mean(xFov,1)),'.'); hold on
hPolFov.MarkerFaceColor = 'k';
hPolFov.MarkerEdgeColor = 'none';
hPolFov.MarkerSize = 1;


[u,v] = pol2cart(angle(mean(xFovAct,1)),log(abs(mean(xFovAct,1))+1));
hPolAct= plot(u,v,'o'); hold on
% hPolAct = polar(angle(mean(xFovAct,1)),abs(mean(xFovAct,1)),'.'); hold on
% hPolAct = polarplot(angle(mean(xFovAct,1)),abs(mean(xFovAct,1)),'.'); hold on
hPolAct.MarkerFaceColor = [1 1 1]-eps;
hPolAct.MarkerEdgeColor = 'none';
hPolAct.MarkerSize = 1;


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



%% Select a single voxel roughly matching the observed group-level effect (just to not confused the reader)
indIn = featSel1.indIn;
act = sqrt(sum(featSel1.featSeq.featVal(:,featSelStep_ind,condPair_ind).^2,2));
actRankScore = abs(prctile(act(featSel1.indIn),100)-act);

xGrat1 = permute(mean(x(y==1,:),1),[2 1]);
xGrat2 = permute(mean(x(y==2,:),1),[2 1]);
xPlaid = permute(mean(x(y==3,:),1),[2 1]);
allDiff   = wrapToPi(angle(mean(cat(2,xGrat1,xGrat2,xPlaid),2)) - angle(mean(cat(1,xGrat1,xGrat2,xPlaid),1)));
gratDiff  =      abs(                 xGrat1                    -                         xGrat2            );
plaidDiff = wrapToPi(angle(                         xPlaid    ) - angle(mean(cat(2,xGrat1,xGrat2       ),2)));
plaidDiffMag =         abs(                         xPlaid    ) -   abs(mean(cat(2,xGrat1,xGrat2       ),2));

indIn = indIn & allDiff./pi*6>-2 & allDiff./pi*6<2;
gratDiffRankScore  = gratDiff;
plaidDiffRankScore = abs(wrapToPi((-0.150/6*pi)-plaidDiff));
plaidDiffMagScore  = abs(plaidDiffMag);

[~,bTmp] = sort(sqrt(...
(0.10.*actRankScore./std(actRankScore      )).^2 + ...
(1.00.*gratDiffRankScore./std(gratDiffRankScore )).^2 + ...
(10.0.*plaidDiffRankScore./std(plaidDiffRankScore)).^2 + ...
(1.00.*plaidDiffMagScore./std(plaidDiffMagScore)).^2 ...
),'ascend');
bTmp = bTmp(indIn);
b = bTmp(1);
voxInd = false(size(indIn));
voxInd(b) = true;

% figure();
% yList = unique(y);
% hPol1vox = {};
% for yInd = 1:length(yList)
%     tmpX = squeeze(x(yList(yInd)==y,b));
%     tmpXmean = mean(tmpX,1);
%     [u,v] = pol2cart(angle(tmpXmean),log(abs(tmpXmean)+1));
%     hPol1vox{yInd} = plot(u,v,'.'); hold on
% end
% for yInd = 1:length(yList)
%     tmpX = squeeze(x(yList(yInd)==y,b));
%     n = size(tmpX,1);
%     tmpXboot = nan(nBoot,size(tmpX,2));
%     for bootInd = 1:nBoot
%         boot = nan(n,1);
%         for i = 1:n
%             boot(i) = randperm(n,1);
%         end
%         tmpXboot(bootInd,:) = mean(tmpX(boot,:),1);
%     end
%     polyCont = credibleInt2D([real(tmpXboot) imag(tmpXboot)],0.05);
%     [theta,rho] = cart2pol(polyCont.Vertices(:,1),polyCont.Vertices(:,2));
%     [polyCont.Vertices(:,1),polyCont.Vertices(:,2)] = pol2cart(theta,log(rho+1));
%     hPol2vox{yInd} = plot(polyCont);
%     hPol2vox{yInd}.LineStyle = 'none';
%     hPol2vox{yInd}.FaceColor = hPol1vox{yInd}.Color;
% end
% axis tight square
% axis(max(abs(axis))*[-1 1 -1 1])
% title(num2str(ii))



%% Plot the selected voxel
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

% legend([hPolFov hPolAct hPol1vox],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid'},'box','off')
legend([hPolFov hPolAct hPol1vox],{'v1 roi vox' 'selected vox' 'grat1' 'grat2' 'plaid'},'box','off')




function [f,voxInd] = plotHr(d,p,featSel1,featSel2,voxInd,visibilityFlag)
nBoot = p.boot.n;
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Get data
[xSin,~,~] = getXYK(d,p,0);
[x,y,~] = getXYK(d,p,1);
x = cat(3,x,x(:,:,1));


xFov       = x(:,featSel2.indIn                  ,:);
xFovNotAct = x(:,featSel2.indIn & ~featSel1.indIn,:);
xFovAct    = x(:,featSel1.indIn                  ,:);
t = (1:size(xFovAct,3))-1;
tFull = linspace(0,13,13*10);

if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

tmpX = squeeze(mean(xFovNotAct,1));
% tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
% tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
hPlotFov = plot(t,tmpX,'Color','k'); hold on
set(hPlotFov,'Color',[get(hPlotFov(1),'Color') 0.15])
set(hPlotFov,'LineWidth',eps);

tmpX = squeeze(mean(xFovAct,1));
% tmpX(tmpX>=0) = log(tmpX(tmpX>=0)+1);
% tmpX(tmpX<0) = -log(-tmpX(tmpX<0)+1);
hPlotAct = plot(t,tmpX,'w'); hold on
set(hPlotAct,'Color',[get(hPlotAct(1),'Color') 0.15])
hPolAct.LineWidth = eps;

ax = gca;
ax.Color = [1 1 1].*0.6;
axis tight


%% Choose one voxel
% featSelSteps_labelList = featSel1.featSeq.featSelList;
% for i = 1:length(featSelSteps_labelList)
%     tmp = strsplit(featSelSteps_labelList{i},':');
%     featSelSteps_labelList(i) = tmp(1);
% end

% featSelStep_ind = ismember(featSelSteps_labelList,{'act' 'respVecSig'});
% condPair_ind = logical([1 0 0 0]);
% if ~all(featSel1.featSeq.condPairList{condPair_ind}==[1 2 3])
%     error('lazy coding')
% end



% act = sqrt(sum(featSel1.featSeq.featVal(:,featSelStep_ind,condPair_ind).^2,2));
% [~,b] = min(abs(prctile(act(featSel1.indIn),80)-act));
% voxInd = false(size(act));
% voxInd(b) = true;


% %% Plot the one voxel
% tmpX = squeeze(x(:,b,:));
% n = size(tmpX,1);
% tmpXmean = mean(tmpX,1);
% tmpXboot = nan(nBoot,size(tmpX,2));
% for bootInd = 1:nBoot
%     boot = nan(n,1);
%     for i = 1:n
%         boot(i) = randperm(n,1);
%     end
%     tmpXboot(bootInd,:) = mean(tmpX(boot,:),1);
% end
% tmpXer = prctile(tmpXboot,[97.5 2.5],1);
% tmpXer(1,:) = tmpXer(1,:) - tmpXmean;
% tmpXer(2,:) = tmpXmean - tmpXer(2,:);
% hPlot2vox = errorbar(t,tmpXmean,tmpXmean-tmpXer(2,:),tmpXer(1,:)-tmpXmean,'ok');
% hPlot2vox.CapSize = 0;
% hPlot2vox.MarkerEdgeColor = 'k';
% hPlot2vox.MarkerFaceColor = 'w';%hPlot2vox.MarkerEdgeColor;



%% Plot conditions for the one voxel
yList = unique(y);
hPlot1vox = {};
cList = colororder;
for yInd = 1:length(yList)
    tmpX = squeeze(x(yList(yInd)==y,voxInd,:));
    tmpXsin = squeeze(xSin(yList(yInd)==y,voxInd));
    tmpXsin = real(tmpXsin)*sin(linspace(0,2*pi,length(tFull))) ...
        + imag(tmpXsin)*cos(linspace(0,2*pi,length(tFull)));

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
    hPlot1vox{yInd} = shadedErrorBar(tFull,tmpXsinMean,tmpXsinEr,'lineprops',{'Color' cList(yInd,:)}); hold on
    delete(hPlot1vox{yInd}.edge)
%     hPlot2vox{yInd} = plot(t,tmpXmean,'Color',cList(yInd,:)); hold on
%     hPlot1vox{yInd}
%     hPlot1vox{yInd} = plot(t,tmpX,'-','Color',cList(yInd,:)); hold on
end
hPlot1vox = [hPlot1vox{:}];

%% FOV and selected FOV voxels
% runTs = squeeze(mean(xFov,2));
% runTs_av = mean(runTs,1);
% runTs_er = bootci(p.boot.n,{@(x)mean(x,1),runTs},'Type','percentile');
% runTs_erNeg = runTs_av - runTs_er(1,:);
% runTs_erPos = runTs_er(2,:) - runTs_av;

% hErr = errorbar(t,runTs_av,runTs_erNeg,runTs_erPos); hold on
% hErr.Marker='o';
% hErr.LineStyle='none';
% hErr.MarkerEdgeColor='k';
% hErr.MarkerFaceColor='k';
% hErr.Color='k';
% hErr.CapSize=0;
% hErr.MarkerSize = 2;
% ylim([-11 11])

runTs = squeeze(mean(xFovAct,2));
runTs_av = mean(runTs,1);
runTs_er = bootci(p.boot.n,{@(x)mean(x,1),runTs},'Type','percentile');
runTs_erNeg = runTs_av - runTs_er(1,:);
runTs_erPos = runTs_er(2,:) - runTs_av;

hErrAct = errorbar(t,runTs_av,runTs_erNeg,runTs_erPos); hold on
hErrAct.Marker='o';
hErrAct.LineStyle='none';
hErrAct.MarkerFaceColor='w';
hErrAct.MarkerEdgeColor='k';
hErrAct.Color='k';
hErrAct.CapSize=0;
hErrAct.MarkerSize = 2;
ylim([-11 11])
xlim([0 12])




%% Decorations
tmpX = mean(xFovAct,1);
ylim([-1 1].*max(abs(tmpX(:))))
% ylim([-1/3 1/3].*max(abs(ylim)))

% uistack(hPlot1vox,'top')
% legend([hPlotFov(1) hPlotAct(1) [hPlot1vox.mainLine] [hPlot1vox.patch] hErr],{'fov vox' 'selected fov vox' 'grat1' 'grat2' 'plaid' '+/-95%CI' '+/-95%CI' '+/-95%CI' 'fov mean +/-95%CI'},'box','off','Location','NorthWest')
legend([hPlotFov(1) hPlotAct(1) [hPlot1vox.mainLine] [hPlot1vox.patch] hErrAct],{'all v1 vox' 'selected vox' 'grat1' 'grat2' 'plaid' '+/-95%CI' '+/-95%CI' '+/-95%CI' 'selected vox mean +/-95%CI'},'box','off','Location','NorthWest')
ax = gca;
ax.Box = 'off';

% ax.XAxis.TickValues = 0:1:11;
ax.XAxis.TickValues = 0:3:12;
