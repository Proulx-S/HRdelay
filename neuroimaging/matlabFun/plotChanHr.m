function f = plotChanHr(chan)

chanCondAList = chan.chanCondAList;
chanCondBList = chan.chanCondBList;
stimList = chan.stimCondList;

chan.hrAv = mean(chan.hr,6);
chan.hrEr = std(chan.hr,[],6)./sqrt(size(chan.hr,6));

chan.hrNormAv = mean(chan.hrNorm,6);
chan.hrNormEr = std(chan.hrNorm,[],6)./sqrt(size(chan.hrNorm,6));


for chanCondAind = 1:length(chanCondAList)
    for chanCondBind = 1:length(chanCondBList)
        if chanCondAind~=chanCondBind
            for sectInd = 1:size(chan.hrNormAv,5)
                if chanCondAind<chanCondBind
                    figure('WindowStyle','docked');
                    % Condition-specific Hr
                    condSpec = squeeze(chan.hrNormAv( chanCondAind,chanCondBind,:,:,sectInd))';
                    nonSpec  = squeeze(chan.hrAv(     chanCondBind,chanCondAind,:,:,sectInd))';
                    
                    
                    av       = nonSpec + condSpec;
                    er       = squeeze(chan.hrNormEr( chanCondAind,chanCondBind,:,:,sectInd))';
                    errorbar(av,er,'CapSize',0); hold on
                    
                    ax = gca;
                    ax.ColorOrderIndex = 1;
                    
                    av       = condSpec;
                    errorbar(av,er,'CapSize',0,'linestyle','--'); hold on
                    
                    
                    
                    av = squeeze(chan.hrNormAv(1,2,:,:,sectInd))';
                    er = squeeze(chan.hrNormEr(1,2,:,:,sectInd))';
                end
            end
        end
    end
end



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


