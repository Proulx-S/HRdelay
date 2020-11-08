function [out,outFixDelay]  = crossCorrOnTrial(ts,tsBase,ind,p,constantFromRun)
plotIt = 0;
for i = 1:size(ind,1)%ha (max R 83; max R negBOLD at delay=1.5 4267; big are and inbetween delay 82)
    curInd = ind(i,:);
    curTs = ts(i,:);
    curtsBase = tsBase(i,:);
    
%     %Exclude first cycle
%     curTs(p.t<p.period) = nan;
%     
    
    
    %Define windows for next step
    p.winStart = (0:p.period:((p.L-1)*p.deltaT))+p.fixDelay;
    p.winStart(1) = []; p.winStart(end) = [];
    if mod(length(p.winStart)/p.nTrialPerDataPoint,2) && length(p.winStart)~=p.nTrialPerDataPoint
        error('number of valid trials is not a multiple of nTrialPerDataPoint')
    end
    p.winStart = p.winStart(1:p.nTrialPerDataPoint:length(p.winStart));
    
    
    for win = 1:length(p.winStart)
        includdeInd = zeros(size(p.t));
        includdeInd((p.t>=p.winStart(win)) & (p.t<p.winStart(win)+p.nTrialPerDataPoint*p.period))=1;
        
        tmpTs = curTs;
        tmpTs(~logical(includdeInd)) = nan;
        
        [corrRes, corrResFixDelay] = runCrossCorr(tmpTs,p,constantFromRun(i),plotIt); % add p.t to plot
        
        %Output best res
        out.constant(i,win) = corrRes.constant;
        out.base(i,win) = curtsBase-out.constant(i,win);
        out.amp(i,win) = corrRes.slope*2/out.base(i,win);
        out.delay(i,win) = p.tShifts(corrRes.lag);
        if corrRes.slope>=0
            out.ampRect(i,win) = out.amp(i,win);
            out.delayRect(i,win) = out.delay(i,win);
        else
            out.ampRect(i,win) = -out.amp(i,win);
            out.delayRect(i,win) = out.delay(i,win)+p.period/2;
        end
        if isfield(corrRes,'RsquaredAdj')
            out.Rsquared(i,win) = corrRes.Rsquared;
        else
            out.Rsquared(i,win) = corrRes.Rsquared;
        end
        out.F(i,win) = corrRes.F;
        out.RMSEall(i,win,:) = corrRes.RMSEall;
        out.RMSE(i,win) = corrRes.RMSE;
        out.RMSEnull(i,win) = corrRes.RMSEnull;
        out.ampRaw(i,win) = corrRes.slope*2;
        out.SNR(i,win) = corrRes.slope*2/corrRes.RMSE;
        out.p(i,win) = corrRes.p;
        
        
        %Output res at fixed delay
        outFixDelay.constant(i,win) = corrResFixDelay.constant;
        outFixDelay.base(i,win) = curtsBase-outFixDelay.constant(i,win);
        outFixDelay.amp(i,win) = corrResFixDelay.slope*2/outFixDelay.base(i,win);
        outFixDelay.delay(i,win) = p.tShifts(corrResFixDelay.lag);
        if isfield(corrResFixDelay,'RsquaredAdj')
            outFixDelay.Rsquared(i,win) = corrResFixDelay.Rsquared;
        else
            outFixDelay.Rsquared(i,win) = corrResFixDelay.Rsquared;
        end
        outFixDelay.F(i,win) = corrResFixDelay.F;
        outFixDelay.RMSE(i,win) = corrResFixDelay.RMSE;
        outFixDelay.ampRaw(i,win) = corrResFixDelay.slope*2;
        outFixDelay.SNR(i,win) = corrResFixDelay.slope*2/corrResFixDelay.RMSE;
        outFixDelay.p(i,win) = corrResFixDelay.p;
    end
end