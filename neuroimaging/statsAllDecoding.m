function statsAllDecoding(p,res,info)


metric = 'auc';
condPairList = info.condPairList;
respFeatList = info.respFeatList;
n = size(res{1,1}.subj.(metric),1);
y = nan(length(condPairList),length(respFeatList),n);
for levelInd = 1:length(respFeatList)
    for condPairInd = 1:length(condPairList)
        y(condPairInd,levelInd,:) = res{condPairInd,levelInd}.subj.(metric);
    end
end

y = permute(y,[2 1 3]);
yInfo = 'respFeat x condPair x subj';
% y = y([2 3 1],:,:);
% respFeatList = respFeatList([2 3 1]);

y = cat(2,y(:,1,:),mean(y(:,2:3,:),2));
% y = cat(1,y(1,:,:),mean(y(2:3,:,:),1));
condPairList{2} = 'gratVSplaid';
condPairList(3) = [];

%% Individual effects
ttest_T = nan(size(y,[1 2]));
ttest_P = nan(size(y,[1 2]));
ttest_FDR1 = nan(size(y,[1 2]));
wlcxn_R = nan(size(y,[1 2]));
wlcxn_P = nan(size(y,[1 2]));
wlcxn_FDR1 = nan(size(y,[1 2]));
ttestAll_T = nan([size(y,1) 1]);
ttestAll_P = nan([size(y,1) 1]);
ttestAll_FDR1 = nan([size(y,1) 1]);
wlcxnAll_R = nan([size(y,1) 1]);
wlcxnAll_P = nan([size(y,1) 1]);
wlcxnAll_FDR1 = nan([size(y,1) 1]);
yAll = sum(y.*[1 2]./3,2);
for levelInd = 1:length(respFeatList)
    [~,ttestAll_P(levelInd,1),~,STATS] = ttest(squeeze(yAll(levelInd,:,:)),0.5,'tail','right');
    ttestAll_T(levelInd,1) = STATS.tstat;
    [wlcxnAll_P(levelInd,1),~,STATS] = signrank(squeeze(yAll(levelInd,:,:)),0.5,'tail','right');
    wlcxnAll_R(levelInd,1) = STATS.signedrank;
    % compare stim cond
    for condPairInd = 1:length(condPairList)
        [~,ttest_P(levelInd,condPairInd),~,STATS] = ttest(squeeze(y(levelInd,condPairInd,:)),0.5,'tail','right');
        ttest_T(levelInd,condPairInd) = STATS.tstat;
        [wlcxn_P(levelInd,condPairInd),~,STATS] = signrank(squeeze(y(levelInd,condPairInd,:)),0.5,'tail','right');
        wlcxn_R(levelInd,condPairInd) = STATS.signedrank;
    end
end

%% Multiple comparisons correction on individual effects
% (1) Two independent response features
ind = ismember(respFeatList,{'cartNoDelay' 'delay'});

tmp = ttest_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttest_FDR1(ind,:) = tmp;
tmp = wlcxn_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxn_FDR1(ind,:) = tmp;

tmp = ttestAll_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttestAll_FDR1(ind,:) = tmp;
tmp = wlcxnAll_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxnAll_FDR1(ind,:) = tmp;

disp('---Two independent response features')
for levelInd = 1:length(respFeatList)
    if ind(levelInd)
        disp(['respFeat[' respFeatList{levelInd} '], condPair[all]'])
        disp(['t test    : T   =' num2str(ttestAll_T(levelInd,1)) ', one-sided p(fdr)=' num2str(ttestAll_P(levelInd,1)) '(' num2str(ttestAll_FDR1(levelInd,1)) ')'])
        disp(['wlcxn test: rank=' num2str(wlcxnAll_R(levelInd,1)) '    , one-sided p(fdr)=' num2str(wlcxnAll_P(levelInd,1)) '(' num2str(wlcxnAll_FDR1(levelInd,1)) ')'])
        disp(' ')
        for condPairInd = 1:length(condPairList)
            disp(['respFeat[' respFeatList{levelInd} '], condPair[' condPairList{condPairInd} ']'])
            disp(['t test    : T   =' num2str(ttest_T(levelInd,condPairInd)) ', one-sided p(fdr)=' num2str(ttest_P(levelInd,condPairInd)) '(' num2str(ttest_FDR1(levelInd,condPairInd)) ')'])
            disp(['wlcxn test: rank=' num2str(wlcxn_R(levelInd,condPairInd)) '    , one-sided p(fdr)=' num2str(wlcxn_P(levelInd,condPairInd)) '(' num2str(wlcxn_FDR1(levelInd,condPairInd)) ')'])
            disp(' ')
        end
    end
end
disp('---')
disp(' ')

% (2) Delay response feature only
ttestAll_FDR2 = nan(size(ttest_P));
wlcxnAll_FDR2 = nan(size(wlcxn_P));
ttest_FDR2 = nan(size(ttest_P));
wlcxn_FDR2 = nan(size(wlcxn_P));
ind = ismember(respFeatList,{'delay'});
tmp = ttest_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
ttest_FDR2(ind,:) = tmp;
tmp = wlcxn_P(ind,:);
tmp(:) = mafdr(tmp(:),'BHFDR','true');
wlcxn_FDR2(ind,:) = tmp;

disp('---(2) Delay response feature only')
for levelInd = 1:length(respFeatList)
    if ind(levelInd)
        disp(['respFeat[' respFeatList{levelInd} '], condPair[all]'])
        disp(['t test    : T   =' num2str(ttestAll_T(levelInd,1)) ', one-sided p(fdr)=' num2str(ttestAll_P(levelInd,1)) '(' num2str(ttestAll_FDR2(levelInd,1)) ')'])
        disp(['wlcxn test: rank=' num2str(wlcxnAll_R(levelInd,1)) '    , one-sided p(fdr)=' num2str(wlcxnAll_P(levelInd,1)) '(' num2str(wlcxnAll_FDR2(levelInd,1)) ')'])
        disp(' ')
        for condPairInd = 1:length(condPairList)
            disp(['respFeat[' respFeatList{levelInd} '], condPair[' condPairList{condPairInd} ']'])
            disp(['t test    : T   =' num2str(ttest_T(levelInd,condPairInd)) ', one-sided p(fdr)=' num2str(ttest_P(levelInd,condPairInd)) '(' num2str(ttest_FDR2(levelInd,condPairInd)) ')'])
            disp(['wlcxn test: rank=' num2str(wlcxn_R(levelInd,condPairInd)) '    , one-sided p(fdr)=' num2str(wlcxn_P(levelInd,condPairInd)) '(' num2str(wlcxn_FDR2(levelInd,condPairInd)) ')'])
            disp(' ')
        end
    end
end
disp('---')
disp(' ')


%%%%%%%%%
%% Anovas
% set intercept to chance
y = y-0.5;

%% ANOVA: Two independent response features
disp('---ANOVA: Two independent response features')
F1 = 'respFeat';
F1L1 = 'cartNoDelay'; F1L1ind = ismember(respFeatList,F1L1);
F1L2 = 'delay'; F1L2ind = ismember(respFeatList,F1L2);
F2 = 'condPair';
F2L1 = 'grat1VSgrat2'; F2L1ind = ismember(condPairList,F2L1);
F2L2 = 'gratVSplaid'; F2L2ind = ismember(condPairList,F2L2);
disp(['ranova: respFeat[' F1L1 ',' F1L2 '] X condPair[' F2L1 ',' F2L2 ']'])
% t = table(squeeze(y(1      ,1,:))      ,squeeze(y(1      ,2      ,:)),squeeze(y(2      ,1      ,:)),squeeze(y(2      ,2      ,:)));
t = table(squeeze(y(F1L1ind,F2L1ind,:)),squeeze(y(F1L1ind,F2L2ind,:)),squeeze(y(F1L2ind,F2L1ind,:)),squeeze(y(F1L2ind,F2L2ind,:)));
t.Properties.VariableNames = {[F1L1 '_' F2L1]...
    [F1L1 '_' F2L2]...
    [F1L2 '_' F2L1]...
    [F1L2 '_' F2L2]};
withinDesign = table(categorical({F1L1 F1L1 F1L2 F1L2}')...
    ,categorical({F2L1 F2L2 F2L1 F2L2}'));
withinDesign.Properties.VariableNames = {F1 F2};
withinModel = [F1 '*' F2];
rm = fitrm(t,[t.Properties.VariableNames{1} '-' t.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign);
tbl = ranova(rm,'WithinModel',withinModel);
disp(tbl(:,[2 4 5]))

% Post-hoc
disp(' ')
disp('post-hoc')
levelList2 = cellstr(unique(withinDesign.respFeat));
tbl = cell(size(levelList2));
for levelInd = 1:length(levelList2)
    ind = ismember(withinDesign.respFeat,levelList2{levelInd});
    t2 = t(:,ind);
    withinDesign2 = withinDesign(ind,2);
    withinModel2 = 'condPair';
    rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
    tbl{levelInd} = ranova(rm,'WithinModel',withinModel2);
    tbl_P(:,levelInd) = tbl{levelInd}(:,5);
    tbl_P.Properties.VariableNames(levelInd) = cellstr(levelList2{levelInd});
end
FDR = nan(size(tbl{levelInd},1),size(tbl_P,2));
for i = 1:2:size(tbl_P,1)
    FDR(i,:) = mafdr(table2array(tbl_P(i,:)),'BHFDR','true');
end
for levelInd = 1:length(levelList2)
    tbl{levelInd} = [tbl{levelInd} table(FDR(:,levelInd),'VariableNames',{'fdr'})];
end
for levelInd = 1:length(levelList2)
    disp(['ranova (respFeat[' levelList2{levelInd} ']): condPair[' F2L1 ',' F2L2 ']'])
    disp(tbl{levelInd}(:,[2 4 5 end]))
    disp(' ')
end

%% ANOVA: Addition of delay
disp('---ANOVA: Addition of delay')
F1 = 'respFeat';
F1L1 = 'cartNoDelay'; F1L1ind = ismember(respFeatList,F1L1);
F1L2 = 'cart'; F1L2ind = ismember(respFeatList,F1L2);
F2 = 'condPair';
F2L1 = 'grat1VSgrat2'; F2L1ind = ismember(condPairList,F2L1);
F2L2 = 'gratVSplaid'; F2L2ind = ismember(condPairList,F2L2);
disp(['ranova: respFeat[' F1L1 ',' F1L2 '] X condPair[' F2L1 ',' F2L2 ']'])
% t = table(squeeze(y(1      ,1,:))      ,squeeze(y(1      ,2      ,:)),squeeze(y(2      ,1      ,:)),squeeze(y(2      ,2      ,:)));
t = table(squeeze(y(F1L1ind,F2L1ind,:)),squeeze(y(F1L1ind,F2L2ind,:)),squeeze(y(F1L2ind,F2L1ind,:)),squeeze(y(F1L2ind,F2L2ind,:)));
t.Properties.VariableNames = {[F1L1 '_' F2L1]...
    [F1L1 '_' F2L2]...
    [F1L2 '_' F2L1]...
    [F1L2 '_' F2L2]};
withinDesign = table(categorical({F1L1 F1L1 F1L2 F1L2}')...
    ,categorical({F2L1 F2L2 F2L1 F2L2}'));
withinDesign.Properties.VariableNames = {F1 F2};
withinModel = [F1 '*' F2];
rm = fitrm(t,[t.Properties.VariableNames{1} '-' t.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign);
tbl = ranova(rm,'WithinModel',withinModel);
disp(tbl(:,[2 4 5]))

% Planned comparison
disp(' ')
disp(['planned comparison: ' F2L2 '; ' F1L1 ' VS ' F1L2])
indF1 = ismember(withinDesign.(F1),{F1L1 F1L2});
indF2 = ismember(withinDesign.(F2),{F2L2});
t2 = t(:,indF1&indF2);
[~,ttest_P,~,STATS] = ttest(table2array(t2(:,1)),table2array(t2(:,2)));
ttest_T = STATS.tstat;
[wlcxn_P,~,STATS] = signrank(table2array(t2(:,1)),table2array(t2(:,2)));
wlcxn_R = STATS.signedrank;
disp(['t test    : T   =' num2str(ttest_T) ', p=' num2str(ttest_P)])
disp(['wlcxn test: rank=' num2str(wlcxn_R) '    , p=' num2str(wlcxn_P)])
disp(' ')







% 
% 
% factorMult = 'condPair';
% factor = 'respFeat';
% % levelList2 = cellstr(unique(withinDesign.(factor)));
% levelList2 = {'gratVSplaid'};
% ind = ismember(withinDesign.respFeat,levelList2{levelInd});
% t2 = t(:,ind);
% withinDesign2 = withinDesign(ind,2);
% withinModel2 = factor;
% 
% 
% tbl = cell(size(levelList2));
% for levelInd = 1:length(levelList2)
%     ind = ismember(withinDesign.respFeat,levelList2{levelInd});
%     t2 = t(:,ind);
%     withinDesign2 = withinDesign(ind,2);
%     withinModel2 = factor;
%     rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
%     tbl{levelInd} = ranova(rm,'WithinModel',withinModel2);
%     tbl_P(:,levelInd) = tbl{levelInd}(:,5);
%     tbl_P.Properties.VariableNames(levelInd) = cellstr(levelList2{levelInd});
% end
% FDR = nan(size(tbl{levelInd},1),size(tbl_P,2));
% for i = 1:2:size(tbl_P,1)
%     FDR(i,:) = mafdr(table2array(tbl_P(i,:)),'BHFDR','true');
% end
% for levelInd = 1:length(levelList2)
%     tbl{levelInd} = [tbl{levelInd} table(FDR(:,levelInd),'VariableNames',{'fdr'})];
% end
% for levelInd = 1:length(levelList2)
%     disp(['ranova (' factorMult '[' levelList2{levelInd} ']): ' factor '[' F2L1 ',' F2L2 ']'])
%     disp(tbl{levelInd}(:,[2 4 5 end]))
%     disp(' ')
% end
% 
% 
% 
% 
% % (1) Two independent response features
% % (2) Delay response feature only
% 
% tbl_P
% 
% 
% 
% 
% 
% 
% 
% 
% for condPairInd = 1:length(condPairList)
%     disp(['ranova: respFeat[' F1L1 ',' F1L2 '] X condPair[' F2L1 ',' F2L2 ']'])
% 
%     disp(['ranova (condPair[' condPairList{condPairInd} ']): respFeat[' F1L1 ',' F1L2 ']'])
%     ind = ~cellfun('isempty',strfind(t.Properties.VariableNames,condPairList{condPairInd}));
%     t2 = t(:,ind);
%     withinDesign2 = withinDesign(ind,1);
%     withinModel2 = 'respFeat';
%     rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
%     tbl = ranova(rm,'WithinModel',withinModel2);
%     disp(tbl(:,[2 4 5]))
% end
% 
% 
% 
% %%
% disp('ranova: respFeat[cartNoDelay,cart] X condPair[grat1VSgrat2,gratVSplaid]')
% t = table(squeeze(y(1,1,:)),squeeze(y(1,2,:)),squeeze(y(3,1,:)),squeeze(y(3,2,:)));
% t.Properties.VariableNames = {'cartNoDelay_grat1VSgrat2'...
%     'cartNoDelay_gratVSplaid'...
%     'cart_grat1VSgrat2'...
%     'cart_gratVSplaid'};
% withinDesign = table(categorical({'cartNoDelay' 'cartNoDelay' 'cart' 'cart'}')...
%     ,categorical({'grat1VSgrat2' 'gratVSplaid' 'grat1VSgrat2' 'gratVSplaid'}'));
% withinDesign.Properties.VariableNames = {'respFeat' 'condPair'};
% withinModel = 'respFeat*condPair';
% rm = fitrm(t,'cartNoDelay_grat1VSgrat2-cart_gratVSplaid~1','WithinDesign',withinDesign);
% tbl = ranova(rm,'WithinModel',withinModel);
% disp(tbl(:,[2 4 5]))
% 
% %Planned comparison
% for condPairInd = 1:length(condPairList)
%     disp(['ranova (condPair[' condPairList{condPairInd} ']): respFeat[cartNoDelay,cart]'])
%     ind = ~cellfun('isempty',strfind(t.Properties.VariableNames,condPairList{condPairInd}));
%     t2 = t(:,ind);
%     withinDesign2 = withinDesign(ind,1);
%     withinModel2 = 'respFeat';
%     rm = fitrm(t2,[t2.Properties.VariableNames{1} '-' t2.Properties.VariableNames{end} '~1'],'WithinDesign',withinDesign2);
%     tbl = ranova(rm,'WithinModel',withinModel2);
%     disp(tbl(:,[2 4 5]))
% end

