function statsChanHr(p,chan)


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
funPath = fullfile(repoPath,'D-tidy\DecodingHR\channelHr');
outDir  = 'stats';
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


if ~exist(fullfile(funPath,outDir),'dir'); mkdir(fullfile(funPath,outDir)); end
diary(fullfile(funPath,outDir,[mfilename '-' datestr(now,'YYYYmmDD-hh_MM_ss') '.stats']))
disp(datestr(now))

for respFeatInd = 1:length(respFeatList)
    wSectFlag = all(~isnan(chan.hr(1,2,1,1,respFeatInd,:,1)));
    for sectInd = 1:length(wSectList)
        if ~isnan(chan.hr(1,2,1,1,respFeatInd,sectInd,1))
            for chanCondAind = 1:length(chanCondAList)
                for chanCondBind = 1:length(chanCondBList)
                    if chanCondAind~=chanCondBind
                        if chanCondAind<chanCondBind
                            %% Non-specific channel (below diag in chan);
                            disp('++++++++++++++++++++++')
                            disp('++++++++++++++++++++++')
                            if wSectFlag
                                disp([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  non-specific channel;  ' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]]);
                            else
                                disp([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  non-specific channel;  ' respFeatList{respFeatInd}]);
                            end
                            disp('++++++++++++++++++++++')
                            disp('++++++++++++++++++++++')
                            x = squeeze(chan.sin(chanCondBind,chanCondAind,:,:,respFeatInd,sectInd,:))';
                            sinRespStats(x)
                            
                            
                            %% Condition-specific channel (above diag in chan);
                            disp('++++++++++++++++++++++')
                            disp('++++++++++++++++++++++')
                            if wSectFlag
                                disp([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  condition-specific channel;  ' [respFeatList{respFeatInd} ' -> ' wSectList{sectInd}]]);
                            else
                                disp([[chanCondAList{chanCondAind} 'VS' chanCondBList{chanCondBind}] ';  condition-specific channel;  ' respFeatList{respFeatInd}]);
                            end
                            disp('++++++++++++++++++++++')
                            disp('++++++++++++++++++++++')
                            x = squeeze(chan.sin(chanCondAind,chanCondBind,:,:,respFeatInd,sectInd,:))';
                            sinRespStats(x)
                        end
                    end
                end
            end
        end
    end
end

diary off





function sinRespStats(x)
% Center delay
theta = wrapToPi(angle(x)-angle(mean(x(:))));
rho = abs(x);
[u,v] = pol2cart(theta,rho);
x = complex(u,v);

% Stats on 2D response vector (conservative)
disp('---------------')
disp('2D Response Vector')
disp('---------------')
disp('Ori1 vs Ori2 vs Plaid (rmanova)')
TBL = manovaX(x);
disp(TBL(2,:))

disp('Grat vs Plaid:')
X = cat(1,mean(x(:,1:2),2),x(:,3));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp('Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f')]);
disp([' p=' num2str(stats.P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = cat(1,x(:,1),x(:,2));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp(' Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f') '; p=' num2str(stats.P,'%0.2f')]);


disp('---------------')
disp('Amplitude')
disp('---------------')

disp('Ori1 vs Ori2 vs Plaid (ranova)')
t = table(abs(x(:,1)),abs(x(:,2)),abs(x(:,3)));
t.Properties.VariableNames = {'grat1' 'grat2' 'plaid'};
withinDesign = table({'grat1' 'grat2' 'plaid'}');
withinDesign.Properties.VariableNames = {'cond'};
withinModel = 'cond';
rm = fitrm(t,'grat1-plaid~1','WithinDesign',withinDesign,'WithinModel',withinModel);
tbl = ranova(rm);
disp(tbl)

disp('Ori1 vs Ori2 vs Plaid (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(abs(x(:,1:3)),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);

disp('Grat vs Plaid (one-tailed):')
X = abs(mean(x(:,1:2),2)); Y = abs(x(:,3));
[~,Po,~,~] = ttest(X,Y,'tail','right');
[~,P,~,STATS] = ttest(X,Y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; one-sided p=' num2str(Po,'%0.2f')  '; two-sided p=' num2str(P,'%0.2f')]);
[Po,~,~] = signrank(X,Y,'tail','right');
[P,~,STATS] = signrank(X,Y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; one-sided p=' num2str(Po,'%0.2f') '; two-sided p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = abs(x(:,1)); Y = abs(x(:,2));
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
    X = abs(x(:,allPairs(pairInd,1))); Y = abs(x(:,allPairs(pairInd,2)));
    [P(pairInd),~,STATS] = signrank(X,Y);
    sRank(pairInd) = STATS.signedrank;
    disp(['cond ' num2str(allPairs(pairInd,:)) ': signed rank=' num2str(sRank(pairInd)) ', p=' num2str(P(pairInd))])    
end



disp('---------------')
disp('Delay')
disp('---------------')

disp('Ori1 vs Ori2 vs Plaid (ranova)')
t = table(angle(x(:,1)),angle(x(:,2)),angle(x(:,3)));
t.Properties.VariableNames = {'grat1' 'grat2' 'plaid'};
withinDesign = table({'grat1' 'grat2' 'plaid'}');
withinDesign.Properties.VariableNames = {'cond'};
withinModel = 'cond';
rm = fitrm(t,'grat1-plaid~1','WithinDesign',withinDesign,'WithinModel',withinModel);
tbl = ranova(rm);
disp(tbl)

disp('Ori1 vs Ori2 vs Plaid (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(angle(x(:,1:3)),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);

disp('Grat vs Plaid (one-tailed):')
X = angle(mean(x(:,1:2),2)); Y = angle(x(:,3));
[~,Po,~,STATS] = ttest(X,Y,'tail','right');
[~,P,~,STATS] = ttest(X,Y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; one-sided p=' num2str(Po,'%0.2f') '; two-sided p=' num2str(P,'%0.2f')]);
[Po,~,STATS] = signrank(X,Y,'tail','right');
[P,~,STATS] = signrank(X,Y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; one-sided p=' num2str(Po,'%0.2f') '; two-sided p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = angle(x(:,1)); Y = angle(x(:,2));
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
    X = angle(x(:,allPairs(pairInd,1))); Y = angle(x(:,allPairs(pairInd,2)));
    [P(pairInd),~,STATS] = signrank(X,Y);
    sRank(pairInd) = STATS.signedrank;
    disp(['cond ' num2str(allPairs(pairInd,:)) ': signed rank=' num2str(sRank(pairInd)) ', p=' num2str(P(pairInd))])
end
