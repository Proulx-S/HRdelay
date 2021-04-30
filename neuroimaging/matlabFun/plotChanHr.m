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
                            subplot(3,3,[4 7])
                            % Show the condition-specific Hrs
                            nonSpecChan_av  = squeeze(chan.hrAv(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd))';
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
                            
                            subplot(3,3,1)
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
                            
                            
                            %% Non-specific channel (below diag in chan);
                            % Hr sin summary
                            x = squeeze(chan.sin(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd,:))';
                            ftmp = sinRespPlot(p,x);
                            subplot(3,3,3,copyobj(ftmp(1).Children(1), f(end)))
                            ftmp(1).Children(3).Title.String = [];
                            subplot(3,3,6,copyobj(ftmp(1).Children(3), f(end)))
                            subplot(3,3,[2 5],copyobj(ftmp(2).Children, f(end)))
                            subplot(3,3,[8 9],copyobj(ftmp(3).Children, f(end)))
                            close(ftmp)
                            
                            if wSectFlag
                                suptitle([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  non-specific channel;  ' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]]);
                            else
                                suptitle([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  non-specific channel;  ' [respFeatList{respFeatInd}]]);
                            end

                            
                            
                            
                            %% Condition-specific channel (above diag in chan)
                            if p.figOption.verbose>1
                                f = [f figure('WindowStyle','docked','visible','on')];
                            elseif p.figOption.verbose>0
                                f = [f figure('WindowStyle','docked','visible','off')];
                            else
                                disp('p.figOption.verbose=0, not plotting anything')
                                return
                            end
                            subplot(3,3,1)
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
                            
                            subplot(3,3,[4 7])
                            % The normalized condition-specific HRs are
                            % hard to interprete. To better visualize the
                            % effect of stimulus condition on the Hr, let's
                            % just add the normalized condition-specific
                            % HRs a generic Hr derived from the
                            % grand-average Hr from the non-specific
                            % channel.
                            genHr  = mean(nonSpecChan_av,2);
                            av       = genHr + condSpecChan_Nav; % add the condSpec response to the nonSpec response
                            er       = condSpecChan_Ner;
                            errorbar(av,er,'CapSize',0); hold on
                            ylabel({'generic Hr + ' 'normalized condition-specific HRs'})
                            legend(stimList,'location','northwest')
                            
                            %% Condition-specific channel (above diag in chan)
                            % Hr sin summary
                            genX = squeeze(chan.sin(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd,:))';
                            genX = mean(genX(:));
                            x = squeeze(chan.sin(chanCondAind,chanCondBind,:,:,respFeatInd,sectInd,:))';
                            x = x - mean(x,2) + genX;
                            
                            ftmp = sinRespPlot(p,x);
                            subplot(3,3,3,copyobj(ftmp(1).Children(1), f(end)))
                            ftmp(1).Children(3).Title.String = [];
                            subplot(3,3,6,copyobj(ftmp(1).Children(3), f(end)))
                            subplot(3,3,[2 5],copyobj(ftmp(2).Children, f(end)))
                            subplot(3,3,[8 9],copyobj(ftmp(3).Children, f(end)))
                            close(ftmp)
                            
                            if wSectFlag
                                suptitle([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  condition-specific channel;  ' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]]);
                            else
                                suptitle([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  condition-specific channel;  ' [respFeatList{respFeatInd}]]);
                            end
                            drawnow
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




function f = sinRespPlot(p,x)
f = [];
colors = [  0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250];
% Remove subject effect
x = x - mean(x,2) + mean(x(:));
% theta = angle(x) - angle(mean(x,2)) + angle(mean(x(:)));
% rho = abs(x) ./ abs(mean(x,2)) .* abs(mean(x(:)));
% [u,v] = pol2cart(theta,rho); clear theta rho
% x = complex(u,v); clear u v

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
ax.ThetaTickLabel = ax.ThetaTick(1:end)/360*12;
% ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
% ax.ThetaTickLabel(1,:) = '0 ';
ax.ThetaAxis.Label.String = {'delay' '(sec)'};
ax.ThetaAxis.Label.Rotation = 0;
ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.RAxis.Label.String = 'amp (%BOLD)';
ax.RAxis.Label.Rotation = 77;

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
Axis = axis;
axis(Axis);
[u,v] = pol2cart(linspace(-pi,pi,100),abs(mean(x(:))));
plot(u,v,'k')
[~,rho] = cart2pol(Axis(1),Axis(4));
[u,v] = pol2cart(angle(mean(x(:))),linspace(0,rho,100));
plot(u,v,'k');
ax = gca;
set(ax.Children(1:2),'color',ones(1,3).*0.5)
uistack(ax.Children(1:2),'bottom');


% Bars
x2 = abs(x) - abs(mean(x(:)));
f = [f figure('WindowStyle','docked')]; clear hB hErr
for condInd = 1:3
    x2av = mean(x2(:,condInd));
    x2er = std(x2(:,condInd))./sqrt(size(x2,1));
    hB(condInd) = bar(condInd,x2av); hold on
    hB(condInd).FaceColor = colors(condInd,:);
    hErr(condInd) = errorbar(condInd,x2av,x2er); hold on
    hErr(condInd).LineStyle = 'none';
    hErr(condInd).Marker = 'none';
    hErr(condInd).Color = 'k';
end
ax = gca;
ylabel({'Response amplitude' '(%BOLD)'})
ax.XTick = 1:3;
ax.XTickLabel = {'grat1' 'grat2' 'plaid'};
box off
amp = [hB.YData];
disp('***')
disp(['amp (plaid-grat) = ' num2str(amp(3)-mean(amp(1:2)),'%0.3f%%BOLD')])
disp(['amp (grat2-grat1) = ' num2str(amp(2)-amp(1),'%0.3f%%BOLD')])
disp('***')

x2 = angle(x) - angle(mean(x(:)));
x2 = wrapToPi(x2)./pi*6;
f = [f figure('WindowStyle','docked')]; clear hB hErr
for condInd = 1:3
    x2av = mean(x2(:,condInd));
    x2er = std(x2(:,condInd))./sqrt(size(x2,1));
    hB(condInd) = barh(condInd,x2av); hold on
    hB(condInd).FaceColor = colors(condInd,:);
    hErr(condInd) = errorbar(x2av,condInd,[],[],x2er,x2er); hold on
%     hErr(condInd) = errorbar(mean(rad2sec(angle(x(:,condInd))),1),condInd,[],[],std(rad2sec(angle(x(:,condInd))),[],1),std(rad2sec(angle(x(:,condInd))),[],1)); hold on
    hErr(condInd).LineStyle = 'none';
    hErr(condInd).Marker = 'none';
    hErr(condInd).Color = 'k';
end
ax = gca;
xlabel({'Response delay' '(sec)'})
ax.YTick = 1:3;
ax.YTickLabel = {'grat1' 'grat2' 'plaid'};
box off
delay = [hB.YData];
hTxt = text(0,hB(3).XData+hB(3).BarWidth/2,['plaid-grat = ' num2str(delay(3)-mean(delay(1:2)),'%0.3fsec')],'VerticalAlignment','bottom');
disp('***')
disp(['delay (plaid-grat) = ' num2str(delay(3)-mean(delay(1:2)),'%0.3fs')])
disp(['delay (ori2-ori1) = ' num2str(delay(2)-delay(1),'%0.3fs')])
disp('***')





