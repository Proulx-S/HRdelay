function plotChannelHr(res)



condList = unique(res.sess.y{1});
tList = 1:size(res.sess.yHat{1},2);
chanHr  = nan([length(condList) length(condList) length(condList) 1 12 size(res.sess.y,1) size(res.sess.y,2)]);
% chanCondA x chanCondB x stimCond x rep x t x subj x sess
for subjInd = 1:size(res.sess.y,1)
    for sessInd = 1:size(res.sess.y,2)
        repN = length(res.sess.y{subjInd,sessInd})/length(condList);

        y     = res.sess.y{subjInd,sessInd};
        spec  = res.sess.yHat{subjInd,sessInd};
        nSpec = res.sess.yHat_nSpec{subjInd,sessInd};
        
        % hr = chanCondA x chanCondB x stimCond x rep x t
        hr = nan(length(condList),length(condList),length(condList),repN,length(tList));
        % add nSpec and spec together into hr
        for chanCondAind = 1:length(condList)
            for chanCondBind = 1:length(condList)
                for stimCondInd = 1:length(condList)
                    if chanCondAind==chanCondBind
                        hr(chanCondAind,chanCondBind,stimCondInd,:,:) = nSpec(y==condList(stimCondInd),:);
                    else
                        hr(chanCondAind,chanCondBind,stimCondInd,:,:) = spec(y==condList(stimCondInd),:);
                    end
                end
            end
        end
        % remove baseline
        hr = hr - mean(hr,5);
        figure('WindowStyle','docked');
        
        % average reps and output
        chanHr(:,:,:,:,:,subjInd,sessInd) = mean(hr,4);
    end
end
% chanHr = chanCondA x chanCondB x stimCond x rep x t x subj x sess
chanHr = mean(chanHr,7);

% average sess
yHat.spec = mean(yHat.spec,7); % chanCondA x chanCondB x stimCond x rep x t x subj
yHat.nSpec = mean(yHat.nSpec,7); % chanCondA x chanCondB x stimCond x rep x t x subj
% summarize
yHat.spec_av = mean(yHat.spec,6);
yHat.spec_er = std(yHat.spec,[],6)./sqrt(size(yHat.spec,6));
yHat.nSpec_av = mean(yHat.nSpec,6);
yHat.nSpec_er = std(yHat.nSpec,[],6)./sqrt(size(yHat.nSpec,6));



% chanCondA x chanCondB x stimCond x rep x t x subj
figure('WindowStyle','docked');
hEr = [];
y_specAv = squeeze(yHat.spec_av(1,2,:,:))';
y_specEr = squeeze(yHat.spec_er(1,2,:,:))';
hEr = [hEr errorbar(y_specAv,y_specEr)]; hold on
y_nSpecAv = squeeze(yHat.nSpec_av(1,2,:,:))';
y_nSpecEr = squeeze(yHat.nSpec_er(1,2,:,:))';
hEr = [hEr errorbar(y_nSpecAv,y_nSpecEr)]; hold on
hEr(3).Color = 'k';


