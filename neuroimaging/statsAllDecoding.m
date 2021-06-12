function statsAllDecoding(p,res,info)


metric = 'auc';
condPairList = info.condPairList;
respFeatList = info.respFeatList;
n = size(res{1,1}.subj.(metric),1);
y = nan(length(condPairList),length(respFeatList),n);
for respFeatInd = 1:length(respFeatList)
    for condPairInd = 1:length(condPairList)
        y(condPairInd,respFeatInd,:) = res{condPairInd,respFeatInd}.subj.(metric);
    end
end

y = permute(y,[2 1 3]);
yInfo = 'respFeat x condPair x subj';
y = y([2 3 1],:,:);
respFeatList = respFeatList([2 3 1]);

y = cat(2,y(:,1,:),mean(y(:,2:3,:),2));
% y = cat(1,y(1,:,:),mean(y(2:3,:,:),1));
condPairList{2} = 'gratVSplaid';
condPairList(3) = [];

%% All bars
ttest_T = nan(size(y,[1 2]));
ttest_P = nan(size(y,[1 2]));
ttest_FDR = nan(size(y,[1 2]));
wlcxn_R = nan(size(y,[1 2]));
wlcxn_P = nan(size(y,[1 2]));
wlcxn_FDR = nan(size(y,[1 2]));
ttest2_T = nan([size(y,1) 1]);
ttest2_P = nan([size(y,1) 1]);
ttest2_FDR = nan([size(y,1) 1]);
wlcxn2_R = nan([size(y,1) 1]);
wlcxn2_P = nan([size(y,1) 1]);
wlcxn2_FDR = nan([size(y,1) 1]);
for respFeatInd = 1:length(respFeatList)
    [~,ttest2_P(respFeatInd,1),~,STATS] = ttest(squeeze(sum(y(respFeatInd,:,:).*[1 2]./3,2)),0.5,'tail','right');
    ttest2_T(respFeatInd,1) = STATS.tstat;
    [wlcxn2_P(respFeatInd,1),~,STATS] = signrank(squeeze(sum(y(respFeatInd,:,:).*[1 2]./3,2)),0.5,'tail','right');
    wlcxn2_R(respFeatInd,1) = STATS.signedrank;
    for condPairInd = 1:length(condPairList)
        [~,ttest_P(respFeatInd,condPairInd),~,STATS] = ttest(squeeze(y(respFeatInd,condPairInd,:)),0.5,'tail','right');
        ttest_T(respFeatInd,condPairInd) = STATS.tstat;
        [wlcxn_P(respFeatInd,condPairInd),~,STATS] = signrank(squeeze(y(respFeatInd,condPairInd,:)),0.5,'tail','right');
        wlcxn_R(respFeatInd,condPairInd) = STATS.signedrank;
    end
end
tmp = ttest_P([1 2],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttest_FDR([1 2],:) = tmp;
tmp = wlcxn_P([1 2],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxn_FDR([1 2],:) = tmp;
tmp = ttest_P([3],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttest_FDR([3],:) = tmp;
tmp = wlcxn_P([3],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxn_FDR([3],:) = tmp;
tmp = ttest2_P([1 2],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttest2_FDR([1 2],:) = tmp;
tmp = wlcxn2_P([1 2],:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxn2_FDR([1 2],:) = tmp;

disp('---')
for respFeatInd = 1:length(respFeatList)
    disp(['respFeat[' respFeatList{respFeatInd} '], condPair[all]'])
    disp(['t test    : T   =' num2str(ttest2_T(respFeatInd,1)) ', one-sided p(fdr)=' num2str(ttest2_P(respFeatInd,1)) '(' num2str(ttest2_FDR(respFeatInd,1)) ')'])
    disp(['wlcxn test: rank=' num2str(wlcxn2_R(respFeatInd,1)) '    , one-sided p(fdr)=' num2str(wlcxn2_P(respFeatInd,1)) '(' num2str(wlcxn2_FDR(respFeatInd,1)) ')'])
    disp(' ')
    for condPairInd = 1:length(condPairList)
        disp(['respFeat[' respFeatList{respFeatInd} '], condPair[' condPairList{condPairInd} ']'])
        disp(['t test    : T   =' num2str(ttest_T(respFeatInd,condPairInd)) ', one-sided p(fdr)=' num2str(ttest_P(respFeatInd,condPairInd)) '(' num2str(ttest_FDR(respFeatInd,condPairInd)) ')'])
        disp(['wlcxn test: rank=' num2str(wlcxn_R(respFeatInd,condPairInd)) '    , one-sided p(fdr)=' num2str(wlcxn_P(respFeatInd,condPairInd)) '(' num2str(wlcxn_FDR(respFeatInd,condPairInd)) ')'])
        disp(' ')
    end
end
disp('---')
disp(' ')

disp('---')
disp('cartNoAmp')
y2 = (y(:,1,:)+y(:,2,:)*2)/3;
[~,ttest_P2,~,STATS] = ttest(squeeze(y2(2,:,:)),0.5,'tail','right');
disp(['ttest (all): T=' num2str(STATS.tstat) ', p=' num2str(mafdr(ttest_P(2,:),'BHFDR','true'))])
disp(['ttest (grat1VSgrat2 gratVSplaid): ' num2str(mafdr(ttest_P(2,:),'BHFDR','true'))])

ttest_T2 = STATS.tstat;
[wlcxn_P2,~,STATS] = signrank(squeeze(y2(2,:,:)),0.5,'tail','right');
wlcxn_R2 = STATS.signedrank;


% set intercept to chance
y = y-0.5;

disp('ranova: respFeat[cartNoDelay,cartNoAmp] X condPair[grat1VSgrat2,gratVSplaid]')
t = table(squeeze(y(1,1,:)),squeeze(y(1,2,:)),squeeze(y(2,1,:)),squeeze(y(2,2,:)));
t.Properties.VariableNames = {'cartNoDelay_grat1VSgrat2'...
    'cartNoDelay_gratVSplaid'...
    'cartNoAmp_grat1VSgrat2'...
    'cartNoAmp_gratVSplaid'};
withinDesign = table(categorical({'cartNoDelay' 'cartNoDelay' 'cartNoAmp' 'cartNoAmp'}')...
    ,categorical({'grat1VSgrat2' 'gratVSplaid' 'grat1VSgrat2' 'gratVSplaid'}'));
withinDesign.Properties.VariableNames = {'respFeat' 'condPair'};
withinModel = 'respFeat*condPair';
rm = fitrm(t,'cartNoDelay_grat1VSgrat2-cartNoAmp_gratVSplaid~1','WithinDesign',withinDesign);
tbl = ranova(rm,'WithinModel',withinModel);
disp(tbl(:,[2 4 5]))


for respFeatInd = 1:length(respFeatList)
    disp(['ranova (respFeat[' respFeatList{respFeatInd} ']): condPair[grat1VSgrat2,gratVSplaid]'])
    ind = ~cellfun('isempty',strfind(t.Properties.VariableNames,respFeatList{respFeatInd}));
    t2 = t(:,ind);
    withinDesign2 = withinDesign(ind,2);
    withinModel2 = 'condPair';
    rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
    tbl = ranova(rm,'WithinModel',withinModel2);
    disp(tbl(:,[2 4 5]))
end


for condPairInd = 1:length(condPairList)
    disp(['ranova (condPair[' condPairList{condPairInd} ']): respFeat[cartNoDelay,cartNoAmp]'])
    ind = ~cellfun('isempty',strfind(t.Properties.VariableNames,condPairList{condPairInd}));
    t2 = t(:,ind);
    withinDesign2 = withinDesign(ind,1);
    withinModel2 = 'respFeat';
    rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
    tbl = ranova(rm,'WithinModel',withinModel2);
    disp(tbl(:,[2 4 5]))
end



%%
disp('ranova: respFeat[cartNoDelay,cart] X condPair[grat1VSgrat2,gratVSplaid]')
t = table(squeeze(y(1,1,:)),squeeze(y(1,2,:)),squeeze(y(3,1,:)),squeeze(y(3,2,:)));
t.Properties.VariableNames = {'cartNoDelay_grat1VSgrat2'...
    'cartNoDelay_gratVSplaid'...
    'cart_grat1VSgrat2'...
    'cart_gratVSplaid'};
withinDesign = table(categorical({'cartNoDelay' 'cartNoDelay' 'cart' 'cart'}')...
    ,categorical({'grat1VSgrat2' 'gratVSplaid' 'grat1VSgrat2' 'gratVSplaid'}'));
withinDesign.Properties.VariableNames = {'respFeat' 'condPair'};
withinModel = 'respFeat*condPair';
rm = fitrm(t,'cartNoDelay_grat1VSgrat2-cart_gratVSplaid~1','WithinDesign',withinDesign);
tbl = ranova(rm,'WithinModel',withinModel);
disp(tbl(:,[2 4 5]))

%Planned comparison
for condPairInd = 1:length(condPairList)
    disp(['ranova (condPair[' condPairList{condPairInd} ']): respFeat[cartNoDelay,cartNoAmp]'])
    ind = ~cellfun('isempty',strfind(t.Properties.VariableNames,condPairList{condPairInd}));
    t2 = t(:,ind);
    withinDesign2 = withinDesign(ind,1);
    withinModel2 = 'respFeat';
    rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
    tbl = ranova(rm,'WithinModel',withinModel2);
    disp(tbl(:,[2 4 5]))
end

