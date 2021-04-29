function f = plotChanHr(chan)

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
                            f = [f figure('WindowStyle','docked')];
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
                            
                            %% Condition-specific channel (above diag in chan)
                            f = [f figure('WindowStyle','docked')];
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




