function f = plotChanHr(p,chan)


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
funPath = fullfile(repoPath,'D-tidy\DecodingHR\channelHr');
outDir  = 'a';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp


%%

chanCondAList = chan.chanCondAList;
chanCondBList = chan.chanCondBList;
stimList      = chan.stimCondList;
wSectList     = chan.wSectList;
respFeatList  = chan.respFeatList;

chan.hrAv = mean(chan.hr,7);
chan.hrEr = std(chan.hr,[],7)./sqrt(size(chan.hr,7));

chan.hrNormAv = mean(chan.hrNorm,7);
chan.hrNormEr = std(chan.hrNorm,[],7)./sqrt(size(chan.hrNorm,7));

f = [];

for respFeatInd = 1:length(respFeatList)
    wSectFlag = all(~isnan(chan.hr(1,2,1,1,respFeatInd,:,1)));
    for sectInd = 1:length(wSectList)
        if ~isnan(chan.hr(1,2,1,1,respFeatInd,sectInd,1))
            for chanCondAind = 1:length(chanCondAList)
                for chanCondBind = 1:length(chanCondBList)
                    if chanCondAind~=chanCondBind
                        if chanCondAind<chanCondBind
                            %% Non-specific channel (below diag in chan);
                            % Hr time courses
                            if p.figOption.verbose>1
                                f = [f figure('WindowStyle','docked','visible','on')];
                            elseif p.figOption.verbose>0
                                f = [f figure('WindowStyle','docked','visible','off')];
                            else
                                disp('p.figOption.verbose=0, not plotting anything')
                                return
                            end
                            subplot(3,1,[2 3])
                            % Show the condition-specific Hrs
                            nonSpecChan_av  = squeeze(chan.hrAv(    chanCondBind,chanCondAind,:,:,respFeatInd,sectInd))';
                            av = nonSpecChan_av;
                            % using the between-subject error of the
                            % cross-condition differences, by computing the error
                            % on the condition-specific Hrs from which the
                            % non-specific Hr was removed.
                            nonSpecChan_Ner = squeeze(chan.hrNormEr(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd))';
                            er = nonSpecChan_Ner;
                            errorbar(av,er,'CapSize',0); hold on
                            ylabel('condition-specific HRs')
                            legend(stimList,'location','northwest')
                            
                            subplot(3,1,1)
                            % Show the normalized condition-specific HRs.
                            % Cross-condition differences are obscured by
                            % the very large non-specific Hr. Subtracting
                            % it on a subject-by-subject basis
                            % (normalizing) offers a better visualization.
                            nonSpecChan_Nav  = squeeze(chan.hrNormAv(    chanCondBind,chanCondAind,:,:,respFeatInd,sectInd))';
                            av = nonSpecChan_Nav;
                            er = nonSpecChan_Ner;
                            errorbar(av,er,'CapSize',0); hold on
                            ylabel({'normalized' 'condition-specific HRs'})
                            
                            if wSectFlag
                                suptitle({[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] 'non-specific channel' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]});
                            else
                                suptitle({[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] 'non-specific channel' [respFeatList{respFeatInd}]});
                            end
                            
                            %% Non-specific channel (below diag in chan);
                            % Hr sin summary
                            x = squeeze(chan.sin(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd,:))';
                            sinRespPlot(p,x)
                            
                            sin  = squeeze(chan.sin(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd,:))';
                            
                            
                            
                            
                            
                            
                            
                            
                            %% Condition-specific channel (above diag in chan)
                            if p.figOption.verbose>1
                                f = [f figure('WindowStyle','docked','visible','on')];
                            elseif p.figOption.verbose>0
                                f = [f figure('WindowStyle','docked','visible','off')];
                            else
                                disp('p.figOption.verbose=0, not plotting anything')
                                return
                            end
                            subplot(3,1,1)
                            % Show the normalized condition-specific HRs. There is
                            % a significant and spurious non-specific component to
                            % the condition-specific HR, which normalizing
                            % (subtracting the cross-condition average) removes.
                            condSpecChan_Nav = squeeze(chan.hrNormAv(chanCondAind,chanCondBind,:,:,respFeatInd,sectInd))';
                            av       = condSpecChan_Nav;
                            % Use the between-subject error of the normalized
                            % condition-specific HRs (the between-subject error in
                            % the differences between conditions).
                            condSpecChan_Ner = squeeze(chan.hrNormEr(chanCondAind,chanCondBind,:,:,respFeatInd,sectInd))';
                            er       = condSpecChan_Ner;
                            errorbar(av,er,'CapSize',0); hold on
                            ylabel({'normalized' 'condition-specific HRs'})
                            
                            subplot(3,1,[2 3])
                            % The normalized condition-specific HRs are hard to
                            % interprete. To better visualize the effect of
                            % stimulus condition on the Hr, let's just add the
                            % normalized condition-specific HRs to the non-specific
                            % Hr of the non-specific channel (i.e. a single generic
                            % Hr derived from our data).
                            genHr  = mean(nonSpecChan_av,2);
                            av       = genHr + condSpecChan_Nav; % add the condSpec response to the nonSpec response
                            er       = condSpecChan_Ner;
                            errorbar(av,er,'CapSize',0); hold on
                            ylabel({'generic Hr + ' 'normalized condition-specific HRs'})
                            legend(stimList,'location','northwest')
                            
                            if wSectFlag
                                suptitle({[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] 'condition-specific channel' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]});
                            else
                                suptitle({[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] 'condition-specific channel' [respFeatList{respFeatInd}]});
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Save figures
if p.figOption.save
    fullpath = fullfile(funPath,outDir);
    if ~exist(fullpath,'dir'); mkdir(fullpath); end
    fullfilename = fullfile(fullpath,mfilename);
    for i = 1:length(f)
        curF = f(i);
        curF.Color = 'none';
        set(findobj(curF.Children,'type','Axes'),'color','none')
        curFile = [fullfilename '_' num2str(i)];
        curExt = 'svg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curF.Color = 'w';
        curExt = 'fig';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curExt = 'jpg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    end
end




function hrPlot(x,figOption)
% Average sessions
x.hr = mean(x.hr,4);
% Remove random subject effet
[x.hr] = normHr(x);

t = (1:size(x.hr,3))-1;
% offset = mean(x.hr(:));
figure('WindowStyle','docked');
for condInd = 1:size(x.hr,2)
    y = squeeze(x.hr(:,condInd,:));
%     y = y-yBase;
%     h1 = plot(t,y); hold on
%     h2 = plot(t,repmat(yBase,size(t))); hold on
%     for i = 1:length(h1)
%         h2(i).Color = h1(i).Color;
%     end
    yMean = mean(y,1);
    yErr = std(y,[],1)./sqrt(size(y,1));
    hErr(condInd) = errorbar(t,yMean,yErr); hold on
end
set(hErr,'CapSize',0)
xlabel('time from stimulus onset (sec)')
ylabel('%BOLD change')
hLeg = legend({'grat1' 'grat2' 'plaid'});
hLeg.Box = 'off';
grid on



function sinRespPlot(p,x)
f = [];
colors = [  0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250];
% Remove subject effect
theta = wrapToPi( angle(x) - angle(mean(x,2)) + angle(mean(x(:))) );
rho = abs(x) ./ abs(mean(x,2)) .* abs(mean(x(:)));
[u,v] = pol2cart(theta,rho);
x = complex(u,v);

% Scatters
f = [f figure('WindowStyle','docked')];
subplot(2,1,1)
for condInd = 1:3
    h(condInd) = polarplot(angle(x(:,condInd)),abs(x(:,condInd)),'o'); hold on
    h(condInd).Color = colors(condInd,:);
    hM(condInd) = polarplot(angle(mean(x(:,condInd),1)),abs(mean(x(:,condInd),1)),'o'); hold on
    hM(condInd).Color = colors(condInd,:);
    hM(condInd).MarkerFaceColor = colors(condInd,:);
    hM(condInd).MarkerEdgeColor = 'w';
end
set(h,'Marker','.')
set(hM,'MarkerEdgeColor','auto')
set(hM,'MarkerFaceColor','w')
set(hM,'LineWidth',2)
uistack(hM,'top')
hl = legend(char({'grat1' 'grat2' 'plaid'}),'Location','east');
hl.Color = 'none';
hl.Box = 'off';
title('Group Response (subject effect removed)')
ax = gca;
ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
ax.ThetaTickLabel(1,:) = '0 ';
ax.ThetaAxis.Label.String = {'delay' '(sec)'};
ax.ThetaAxis.Label.Rotation = 0;
ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.RAxis.Label.String = 'amp (%BOLD)';
ax.RAxis.Label.Rotation = 80;

subplot(2,1,2)
for condInd = 1:3
    h(condInd) = scatter(real(x(:,condInd)),imag(x(:,condInd)),'.'); hold on
    h(condInd).CData = colors(condInd,:);
    hM(condInd) = scatter(real(mean(x(:,condInd),1)),imag(mean(x(:,condInd),1)),'o'); hold on
    hM(condInd).CData = colors(condInd,:);
    hM(condInd).MarkerFaceColor = 'w';
    hM(condInd).MarkerEdgeColor = colors(condInd,:);
end
set(hM,'LineWidth',2)
ax = gca;
ax.PlotBoxAspectRatio = [1 1 1];
ax.DataAspectRatio = [1 1 1];
axis(axis);
[u,v] = pol2cart(linspace(-pi,pi,100),abs(mean(x(:))));
plot(u,v,'k')
[~,rho] = cart2pol(Axis(1),Axis(4));
[u,v] = pol2cart(angle(mean(x(:))),linspace(0,rho,100));
plot(u,v,'k');


% Bars
f = [f figure('WindowStyle','docked')]; clear hB hErr
for condInd = 1:3
    hB(condInd) = bar(condInd,abs(mean(x(:,condInd),1))); hold on
    hB(condInd).FaceColor = colors(condInd,:);
    hErr(condInd) = errorbar(condInd,abs(mean(x(:,condInd),1)),std(abs(x(:,condInd)),[],1)); hold on
    hErr(condInd).LineStyle = 'none';
    hErr(condInd).Marker = 'none';
    hErr(condInd).Color = 'k';
end
ax = gca;
ylabel({'Response amplitude' '(%BOLD)'})
ax.XTick = 1:3;
ax.XTickLabel = {'grat1' 'grat2' 'plaid'};
f = gcf;
f.Color = 'w';
box off
amp = [hB.YData];
disp('***')
disp(['amp (plaid-grat) = ' num2str(amp(3)-mean(amp(1:2)),'%0.3f%%BOLD')])
disp(['amp (grat2-grat1) = ' num2str(amp(2)-amp(1),'%0.3f%%BOLD')])
disp('***')

f = [f figure('WindowStyle','docked')]; clear hB hErr
for condInd = 1:3
    hB(condInd) = barh(condInd,mean(rad2sec(angle(x(:,condInd))),1)); hold on
    hB(condInd).FaceColor = colors(condInd,:);
    hErr(condInd) = errorbar(mean(rad2sec(angle(x(:,condInd))),1),condInd,[],[],std(rad2sec(angle(x(:,condInd))),[],1),std(rad2sec(angle(x(:,condInd))),[],1)); hold on
    hErr(condInd).LineStyle = 'none';
    hErr(condInd).Marker = 'none';
    hErr(condInd).Color = 'k';
end
ax = gca;
xlabel({'Response delay' '(sec)'})
ax.YTick = 1:3;
ax.YTickLabel = {'grat1' 'grat2' 'plaid'};
f = gcf;
f.Color = 'w';
box off
delay = [hB.YData];
hTxt = text(0,hB(3).XData+hB(3).BarWidth/2,['plaid-grat = ' num2str(delay(3)-mean(delay(1:2)),'%0.3fsec')],'VerticalAlignment','bottom');
disp('***')
disp(['delay (plaid-grat) = ' num2str(delay(3)-mean(delay(1:2)),'%0.3fs')])
disp(['delay (ori2-ori1) = ' num2str(delay(2)-delay(1),'%0.3fs')])
disp('***')


function thetaSec = rad2sec(thetaRad)
thetaSec = -thetaRad/pi*6;

function sinRespStats(x)
% Average sessions
x.sin = mean(x.sin,4);
% Center delay
theta = wrapToPi(angle(x.sin)-angle(mean(x.sin(:))));
rho = abs(x.sin);
[u,v] = pol2cart(theta,rho);
x.sin = complex(u,v);

% Stats on 2D response vector (conservative)
disp('---------------')
disp('2D Response Vector')
disp('---------------')
disp('Grat vs Plaid:')
X = cat(1,mean(x.sin(:,1:2),2),x.sin(:,3));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp('Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f')]);
disp([' p=' num2str(stats.P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = cat(1,x.sin(:,1),x.sin(:,2));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp(' Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f') '; p=' num2str(stats.P,'%0.2f')]);


disp('---------------')
disp('Amplitude')
disp('---------------')

disp('Ori1 vs Ori2 vs Plaid (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(abs(x.sin(:,1:3)),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);

disp('Grat vs Plaid (one-tailed):')
X = abs(mean(x.sin(:,1:2),2)); Y = abs(x.sin(:,3));
[~,P,~,STATS] = ttest(X,Y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; one-sided p=' num2str(P,'%0.2f')]);
[P,~,STATS] = signrank(X,Y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; one-sided p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = abs(x.sin(:,1)); Y = abs(x.sin(:,2));
[~,P,~,STATS] = ttest(X,Y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,~,STATS] = signrank(X,Y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);

allPairs = [1 2; 1 3; 2 3];
P = nan(size(allPairs,1),1);
sRank = nan(size(allPairs,1),1);
for pairInd = 1:size(allPairs,1)
    X = abs(x.sin(:,allPairs(pairInd,1))); Y = abs(x.sin(:,allPairs(pairInd,2)));
    [P(pairInd),~,STATS] = signrank(X,Y);
    sRank(pairInd) = STATS.signedrank;
end



disp('---------------')
disp('Delay')
disp('---------------')

disp('Ori1 vs Ori2 vs Plaid (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(angle(x.sin(:,1:3)),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);

disp('Grat vs Plaid (one-tailed):')
X = angle(mean(x.sin(:,1:2),2)); Y = angle(x.sin(:,3));
[~,P,~,STATS] = ttest(X,Y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; one-sided p=' num2str(P,'%0.2f')]);
[P,~,STATS] = signrank(X,Y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; one-sided p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = angle(x.sin(:,1)); Y = angle(x.sin(:,2));
[~,P,~,STATS] = ttest(X,Y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,~,STATS] = signrank(X,Y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);

allPairs = [1 2; 1 3; 2 3];
P = nan(size(allPairs,1),1);
sRank = nan(size(allPairs,1),1);
for pairInd = 1:size(allPairs,1)
    X = angle(x.sin(:,allPairs(pairInd,1))); Y = angle(x.sin(:,allPairs(pairInd,2)));
    [P(pairInd),~,STATS] = signrank(X,Y);
    sRank(pairInd) = STATS.signedrank;
end


