function YTicks = makeOverlay(axBak,ax,alphaData,cMap,cBar,cLim)

linkaxes([axBak ax]);
linkprop([axBak ax],'position');

hIm = findobj(ax.Children,'Type','Image');
ind = logical(hIm.AlphaData);
hIm.AlphaData(ind) = alphaData.*1;

if exist('cMap','var') && ~isempty(cMap)
    ax.Colormap = cMap;
end
if exist('cLim','var') && ~isempty(cLim)
    ax.CLim = cLim;
else
    cLim = ax.CLim;
end

if exist('cBar','var') && ~isempty(cBar)
    hIm.AlphaData(~ind) = 1;
    ax.YAxis.Visible = 'on';
    
    switch cBar
        case 'log'
            delta = 0.1;
            [YTicks,YTickLabel] = logTicks(ax,delta);
            if length(YTicks)<=2
                delta = 0.01;
                [YTicks,YTickLabel] = logTicks(ax,delta);
                tmp = round(YTicks./delta/10)==YTicks./delta/10;
                YTickLabel = YTicks(tmp);
            end
            
%             YTicks = [];
%             YTickLabel = [];
%             while ~done
%                 YTicks = [YTicks (0:0.1:1)*(10^i)];
%                 YTickLabel = [YTickLabel (10^i)];
%                 i = i+1;
%                 done = YTicks(end)>exp(ax.CLim(2));
%             end
%             imMin = ax.CLim(1);
%             imMax = ax.CLim(2);
%             YTicks = sort(unique(YTicks));
%             YTicks = YTicks(YTicks>exp(imMin) & YTicks<exp(imMax));
            
            
            imMin = ax.CLim(1);
            imMax = ax.CLim(2);
            ax.YTick = interp1([imMin imMax],[1 size(hIm.CData,1)],log(YTicks));
            ax.YTickLabel = cellstr(num2str(YTicks'));
            
%             tmp = ax.YTickLabel;
%             tmp(~tickLabelInd) = {''};
%             ax.YTickLabel = tmp;
            
            tmp = ax.YTickLabel;
            tmp(~ismember(YTicks,YTickLabel)) = {''};
            tmp([1 end]) = cellstr(num2str(str2num(char(ax.YTickLabel([1 end]))),'%0.1f'));
            ax.YTickLabel = tmp;
%             ax.YTickLabel(1) = {num2str(YTicks(1))};
%             ax.YTickLabel(end) = {num2str(YTicks(end))};
            ax.TickDir = 'out';
        case 'lin'
            a = cLim./(abs(cLim) .* 10 .^ (-floor(log10(abs(cLim)))));
            if length(cLim(1):a(1):cLim(2))<10
                a = a(1);
            else
                a = a(2);
            end
            YTicks = 0:a:cLim(2);
            YTicks = YTicks(YTicks>cLim(1) & YTicks<cLim(2));
            if length(YTicks)==2
                YTicks = 0:a/2:cLim(2);
                YTicks = YTicks(YTicks>cLim(1) & YTicks<cLim(2));
            end
            ax.YTick = interp1(cLim,[1 size(hIm.CData,1)],YTicks);
            ax.YTickLabel = cellstr(num2str(YTicks'));
            ax.TickDir = 'out';
        otherwise
            error('X')
    end
end

function [YTicks,YTickLabel] = logTicks(ax,delta)
YTicks = [];
YTickLabel = [];
i = 0;
done = 0;
while ~done
    YTicks = [YTicks (0:delta:1)*(10^i)];
%     YTickLabel = [YTickLabel ((delta*100)^i)];
    YTickLabel = [YTickLabel (10^i)];
    i = i+1;
    done = YTicks(end)>exp(ax.CLim(2));
end
imMin = ax.CLim(1);
imMax = ax.CLim(2);
YTicks = sort(unique(YTicks));
YTicks = YTicks(YTicks>exp(imMin) & YTicks<exp(imMax));
YTicks = [exp(imMin) YTicks exp(imMax)];