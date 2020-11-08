function curResp = interpolateNaverageTimSeries2(splitData_detrended,delay,grandRefDelay,runNum,runList,sess,subjInd,tList,results)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[x,y,z] = ind2sub(size(results.OLS.mixed.Fsess1.val.F),find(results.OLS.mixed.Fsess1.val.F==max(results.OLS.mixed.Fsess1.val.F(:))));
runInd = runList(runNum);
%Shift voxels time series
for x = 1:size(splitData_detrended{runInd},1) %51
    for y = 1:size(splitData_detrended{runInd},2) %15
        for z = 1:size(splitData_detrended{runInd},3) %12
            if any(isnan(delay{subjInd}(x,y,z,sess,:)))
                tmpTs_shift = nan(size(squeeze(splitData_detrended{runInd}(x,y,z,:))));
            else
                tmpTs = squeeze(splitData_detrended{runInd}(x,y,z,:));
                t0 = 0:length(tmpTs)-1;
                tmpTs = [tmpTs(end-11:end); tmpTs; tmpTs(1:12)];
                t0 = [-12:-1 t0 t0(end-11:end)+12];
%                 tShift = wrapToPi(circ_mean(delay{subjInd}(x,y,z,sess,1:2),[],5)+grandRefDelay)/pi*6;
                tShift = wrapToPi(circ_mean(delay{subjInd}(x,y,z,sess,:),[],5)+grandRefDelay)/pi*6;
                %                         figure('WindowStyle','docked');
                %                         plot(t0(13:end-11),tmpTs(13:end-11)); hold on
                
                %Shift and wrap
                t1 = t0+(tShift-round(tShift));
                wrapLength = round(tShift);
                if wrapLength==0
                elseif wrapLength>0
                    t0 = [(-wrapLength:-1)+t0(1) t0];
                    t0(end-(wrapLength-1):end) = [];
                    
                    tmpTsx = tmpTs(end-(wrapLength-1):end);
                    tmpTs = [tmpTsx; tmpTs]; tmpTs(end-(wrapLength-1):end) = [];
                elseif wrapLength<0
                    wrapLength = abs(wrapLength);
                    t0(1:wrapLength) = [];
                    t0(end+1:end+wrapLength) = t0(end)+1:t0(end)+wrapLength;
                    
                    tmpTsx = tmpTs(1:wrapLength);
                    tmpTs = [tmpTs; tmpTsx]; tmpTs(1:wrapLength) = [];
                end
                
                %Interpolate
                tmpTs_shift = interp1(t0,tmpTs,t1,'PCHIP');
                
                %                         t0 = t0+tShift;
                %                         plot(t0(13:end-11),tmpTs_shift(13:end-11)); legend({'orig' 'shift'})
                %                         xlim([-25 125]); ylim([-100 150])
                
                %Cut the extra bit
                t0 = t0(13:end-12);
                tmpTs_shift = tmpTs_shift(13:end-12);
                
                %Wrap back
                wrapLength = -t0(1);
                if wrapLength==0
                elseif wrapLength>0
                    t0 = t0+wrapLength;
                    tmpTs_shiftx = tmpTs_shift(1:wrapLength);
                    tmpTs_shift = [tmpTs_shift tmpTs_shiftx]; tmpTs_shift(1:wrapLength) = [];
                elseif wrapLength<0
                    wrapLength = abs(wrapLength);
                    t0 = t0-wrapLength;
                    tmpTs_shiftx = tmpTs_shift(end-(wrapLength-1):end);
                    tmpTs_shift = [tmpTs_shiftx tmpTs_shift]; tmpTs_shift(end-(wrapLength-1):end) = [];
                end
                
                %                         t0
                %
                %                                                 figure('WindowStyle','docked');
                %                                                 plot(t0,tmpTs_shift)
                %                                                 xlim([-25 125]); ylim([-100 150])
            end
            splitData_detrended{runInd}(x,y,z,:) = tmpTs_shift;
        end
    end
end


curResp = zeros(size(splitData_detrended{1},1),size(splitData_detrended{1},2),size(splitData_detrended{1},3),12); % X x Y x Z x time
for i = 1:length(tList)
    curResp = curResp+splitData_detrended{runInd}(:,:,:,tList(i):(tList(i)+11));
end
curResp = curResp./length(tList);


end

