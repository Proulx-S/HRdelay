function plotChannelHr(res)


%% Compile Hr
condList = unique(res.sess.y{1});
tList = (1:size(res.sess.yHat{1},2))-1;
t = tList'./12.*2.*pi;
X = [sin(t) cos(t) ones(size(t))];
chan.hr     = nan([length(condList) length(condList) length(condList) length(tList) size(res.sess.y,1) size(res.sess.y,2)]);
chan.sin    = nan([length(condList) length(condList) length(condList) 1             size(res.sess.y,1) size(res.sess.y,2)]);
chan.base   = nan([length(condList) length(condList) length(condList) 1             size(res.sess.y,1) size(res.sess.y,2)]);
chan.info1  = 'chanCondA x chanCondB x stimCond x t x subj x sess';
chan.info2  = 'chanCondA x chanCondB: above diag -> cond-specific channel; below diag -> non-specific channel';
chan.info3  = 'cond([grat1VSgrat2 grat1VSplaid grat2VSplaid])';
for subjInd = 1:size(res.sess.y,1)
    for sessInd = 1:size(res.sess.y,2)
        
        y     = res.sess.y{subjInd,sessInd};
        spec  = res.sess.yHat{subjInd,sessInd};
        nSpec = res.sess.yHat_nSpec{subjInd,sessInd};
        
        % hr = chanCondA x chanCondB x stimCond x t
        for chanCondAind = 1:length(condList)
            for chanCondBind = 1:length(condList)
                for stimCondInd = 1:length(condList)
                    if chanCondAind<chanCondBind % upper triangle -> spec
                        Y = mean(spec(y==condList(stimCondInd),:),1);
                    elseif chanCondAind>chanCondBind % lower triangle -> nSpec
                        Y = mean(nSpec(y==condList(stimCondInd),:),1);
                    end
                    % fit sin
                    if chanCondAind~=chanCondBind
                        beta = Y/X';
                        % put in complex format
                        chan.sin(chanCondAind,chanCondBind,stimCondInd,1,subjInd,sessInd) = complex(beta(1),beta(2));
                        % remove baseline and output
                        chan.hr(chanCondAind,chanCondBind,stimCondInd,:,subjInd,sessInd) = Y - beta(3);
                    end
                end
            end
        end
    end
end
%% Average sessions
chan.hr   = mean(chan.hr,6);
chan.sin  = mean(chan.sin,6);

%% Remove subject effect on delay and amplitude
for chanCondAind = 1:length(condList)
    for chanCondBind = 1:length(condList)
        if chanCondAind~=chanCondBind
            hr = permute(chan.hr(chanCondAind,chanCondBind,:,:,:),[4 5 3 1 2]);
            if chanCondAind<chanCondBind
                % when processing cond-specific hr (upper triangle),
                % normalize according to the non-specific hr (lower triangle)
                hrSin = permute(chan.sin(chanCondBind,chanCondAind,:,:,:),[4 5 3 1 2]);
            elseif chanCondAind>chanCondBind
                % when processing non-specific hr (lower triangle),
                % normalize according to itself
                hrSin = permute(chan.sin(chanCondAind,chanCondBind,:,:,:),[4 5 3 1 2]);
            end
            hr = normHr(hr,hrSin);
            chan.hr(chanCondAind,chanCondBind,:,:,:) = permute(hr,[4 5 3 1 2]);
        end
    end
end

%% Remove non-specific response
chan.hrNorm = chan.hr - mean(chan.hr,3);

%% Summarize
chan.hrAv = mean(chan.hr,5);
chan.hrEr = std(chan.hr,[],5);

chan.hrNormAv = mean(chan.hrNorm,5);
chan.hrNormEr = std(chan.hrNorm,[],5);


%% Plot
figure('WindowStyle','docked');
av = squeeze(chan.hrNormAv(1,2,:,:))';
er = squeeze(chan.hrNormEr(1,2,:,:))';
errorbar(av,er); hold on
ax = gca;
ax.ColorOrderIndex = 1;
av = squeeze(chan.hrAv(2,1,:,:))';
er = squeeze(chan.hrEr(2,1,:,:))';
errorbar(av,er,'linestyle','--')


