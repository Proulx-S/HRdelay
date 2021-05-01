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

