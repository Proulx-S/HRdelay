% function xData = finalPrep(data,keepInfo)
% switch keepInfo
%     case {'m','r','i'}
%         complexFlag = 0;
%         xData.data = data.amp;
%         xData.label = data.label;
%     case 'p'
%         complexFlag = 1;
%         [X,Y] = pol2cart(data.delay,1);
%         xData.data = [X Y];
%         xData.label = data.label;
%     case 'mp'
%         complexFlag = 1;
%         [X,Y] = pol2cart(data.delay,data.amp);
%         xData.data = [X Y];
%         xData.label = data.label;
% end
