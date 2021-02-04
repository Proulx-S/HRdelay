function compareRes(X,Y)

metric = 'acc';
metric_lowCI = 'acc_CI5';
metric_highCI = 'acc_CI95';
x = X.sess.(metric)';
y = Y.sess.(metric)';
yNeg = y - Y.sess.(metric_lowCI)';
yPox = Y.sess.(metric_highCI)' - y;
xNeg = x - X.sess.(metric_lowCI)';
xPos = X.sess.(metric_highCI)' - x;
hEb = errorbar(x,y,yNeg,yPox,xNeg,xPos,'o'); hold on
for i = 1:length(hEb)
    hEb(i).MarkerFaceColor = hEb(i).Color;
    hEb(i).CapSize = 0;
    hEb(i).MarkerEdgeColor = 'k';
end

x = X.group.(metric)';
y = Y.group.(metric)';
yNeg = y - Y.group.(metric_lowCI)';
yPox = Y.group.(metric_highCI)' - y;
xNeg = x - X.group.(metric_lowCI)';
xPos = X.group.(metric_highCI)' - x;
hEbM = errorbar(x,y,yNeg,yPox,xNeg,xPos,'^'); hold on
hEbM.MarkerFaceColor = 'w';
hEbM.Color = 'k';

switch metric
    case {'acc' 'auc'}
        xlim([0 1])
        ylim([0 1])
        hChance(1) = plot(xlim,[1 1].*0.5,'-k');
        hChance(2) = plot([1 1].*0.5,ylim,'-k');
        hChance(3) = plot([0 1],[0 1],'-k');
        uistack(hChance,'bottom')
        ax = gca;
        ax.PlotBoxAspectRatio = [1 1 1];
    case 'distT'
        plot(xlim,[1 1].*0,'-k')
        plot([1 1].*0,ylim,'-k')
end
