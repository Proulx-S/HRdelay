function [out,outFixDelay]  = crossCorrOnRun(ts,tsBase,ind,p)
plotIt = 0;
for i = 1:size(ind,1)%ha (max R 83; max R negBOLD at delay=1.5 4267; big are and inbetween delay 82)
    curInd = ind(i,:);
    curTs = ts(i,:);
    curtsBase = tsBase(i,:);
    
    %Exclude first cycle
    curTs(p.t<p.period) = nan;
    
    %Run crossCorr
    [corrRes, corrResFixDelay] = runCrossCorr(curTs,p,[],plotIt); % add p.t to plot
    
    %Output best res
    out.constant(i) = corrRes.constant;
    out.base(i) = curtsBase-out.constant(i);
    out.amp(i) = corrRes.slope*2/out.base(i);
    out.delay(i) = p.tShifts(corrRes.lag);
    if corrRes.slope>=0
        out.ampRect(i) = out.amp(i);
        out.delayRect(i) = out.delay(i);
    else
        out.ampRect(i) = -out.amp(i);
        out.delayRect(i) = out.delay(i)+p.period/2;
    end
    if isfield(corrRes,'RsquaredAdj')
        out.Rsquared(i) = corrRes.Rsquared;
    else
        out.Rsquared(i) = corrRes.Rsquared;
    end
    out.F(i) = corrRes.F;
    out.RMSEall(i,:) = corrRes.RMSEall;
    out.RMSE(i) = corrRes.RMSE;
    out.RMSEnull(i) = corrRes.RMSEnull;
    out.ampRaw(i) = corrRes.slope*2;
    out.SNR(i) = corrRes.slope*2/corrRes.RMSE;
    out.p(i) = corrRes.p;
    
    
    %Output res at fixed delay
    outFixDelay.constant(i) = corrResFixDelay.constant;
    outFixDelay.base(i) = curtsBase-outFixDelay.constant(i);
    outFixDelay.amp(i) = corrResFixDelay.slope*2/outFixDelay.base(i);
    outFixDelay.delay(i) = p.tShifts(corrResFixDelay.lag);
    if isfield(corrResFixDelay,'RsquaredAdj')
        outFixDelay.Rsquared(i) = corrResFixDelay.Rsquared;
    else
        outFixDelay.Rsquared(i) = corrResFixDelay.Rsquared;
    end
    outFixDelay.F(i) = corrResFixDelay.F;
    outFixDelay.RMSE(i) = corrResFixDelay.RMSE;
    outFixDelay.ampRaw(i) = corrResFixDelay.slope*2;
    outFixDelay.SNR(i) = corrResFixDelay.slope*2/corrResFixDelay.RMSE;
    outFixDelay.p(i) = corrResFixDelay.p;    
end