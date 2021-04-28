function f = plotChannelHr(res)


%% Compile Hr
condList = unique(res.sess.y{1});
tList = (1:size(res.sess.yHat{1},2))-1;
t = tList'./12.*2.*pi;
X = [sin(t) cos(t) ones(size(t))];
chan.hr     = nan([length(condList) length(condList) length(condList) length(tList) size(res.sess.y,1) size(res.sess.y,2) size(res.sess.yHat{1},3)]);
chan.sin    = nan([length(condList) length(condList) length(condList) 1             size(res.sess.y,1) size(res.sess.y,2) size(res.sess.yHat{1},3)]);
chan.wSectList = permute({'+real' '+real+imag' '+imag' '-real+imag'},[1 3 4 5 6 7 2]);
chan.info1  = 'chanCondA x chanCondB x stimCond x t x subj x sess x wSect';
chan.info2  = 'chanCondA x chanCondB: above diag -> cond-specific channel; below diag -> non-specific channel';
chan.info3  = 'cond([grat1VSgrat2 grat1VSplaid grat2VSplaid])';
for subjInd = 1:size(res.sess.y,1)
    for sessInd = 1:size(res.sess.y,2)
        
        y     = res.sess.y{subjInd,sessInd};
        for sectInd = 1:size(res.sess.yHat{subjInd,sessInd},3)
            spec  = res.sess.yHat{subjInd,sessInd}(:,:,sectInd);
            nSpec = res.sess.yHat_nSpec{subjInd,sessInd}(:,:,sectInd);
            
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
                            % put sin fit in complex format
                            chan.sin(chanCondAind,chanCondBind,stimCondInd,1,subjInd,sessInd,sectInd) = complex(beta(1),beta(2));
                            % remove baseline from hr
                            chan.hr( chanCondAind,chanCondBind,stimCondInd,:,subjInd,sessInd,sectInd) = Y - beta(3);
                        end
                    end
                end
            end
        end
    end
end
%% Average sessions
chan.hr   = mean(chan.hr,6);
chan.sin  = mean(chan.sin,6);

% Average sin over wSect since we will never use them individually
chan.sin = mean(chan.sin,7);

%% Remove subject effect on delay and amplitude
for chanCondAind = 1:length(condList)
    for chanCondBind = 1:length(condList)
        if chanCondAind~=chanCondBind
            for sectInd = 1:size(chan.hr,7)
                hr = permute(chan.hr(chanCondAind,chanCondBind,:,:,:,:,sectInd),[4 5 3 1 2]);
                if chanCondAind<chanCondBind
                    % when processing cond-specific hr (upper triangle),
                    % normalize according to the non-specific hr (lower triangle)
                    hrSin = permute(chan.sin(chanCondBind,chanCondAind,:,:,:,:),[4 5 3 1 2]);
                elseif chanCondAind>chanCondBind
                    % when processing non-specific hr (lower triangle),
                    % normalize according to itself
                    hrSin = permute(chan.sin(chanCondAind,chanCondBind,:,:,:,:),[4 5 3 1 2]);
                end
                hr = normHr(hr,hrSin);
                chan.hr(chanCondAind,chanCondBind,:,:,:,:,sectInd) = permute(hr,[4 5 3 1 2]);
            end
        end
    end
end

%% Remove non-specific response
chan.hrNorm = chan.hr - mean(chan.hr,3);

%% Summarize
chan.hr     = permute(chan.hr,[1 2 3 4 7 5 6]);
chan.hrNorm = permute(chan.hrNorm,[1 2 3 4 7 5 6]);
chan.sin    = permute(chan.sin,[1 2 3 4 7 5 6]);
chan.wSectList    = permute(chan.wSectList,[1 2 3 4 7 5 6]);
chan.info1  = 'chanCondA x chanCondB x stimCond x t x wSect x subj x sess';


chan.hrAv = mean(chan.hr,6);
chan.hrEr = std(chan.hr,[],6)./sqrt(size(chan.hr,6));

chan.hrNormAv = mean(chan.hrNorm,6);
chan.hrNormEr = std(chan.hrNorm,[],6)./sqrt(size(chan.hrNorm,6));


%% Plot
f = cell(size(chan.hrNormAv,5));
for sectInd = 1:size(chan.hrNormAv,5)
    f{sectInd} = figure('WindowStyle','docked');
    av = squeeze(chan.hrNormAv(1,2,:,:,sectInd))';
    er = squeeze(chan.hrNormEr(1,2,:,:,sectInd))';
    errorbar(av,er,'CapSize',0); hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    av = squeeze(chan.hrAv(2,1,:,:,sectInd))';
    er = squeeze(chan.hrEr(2,1,:,:,sectInd))';
    errorbar(av,er,'linestyle','--','CapSize',0)
    title(chan.wSectList{sectInd})
end


