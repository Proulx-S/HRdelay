function f = showGLMdesign(GLMs)

f = {};


f = [f {figure('WindowStyle','docked')}];
curDesign = GLMs.random.sinDesign;
imagesc(curDesign,[-1 1]); colormap gray
title('Sinusoidal Response');
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
ylabel('TRs');
ax1.XTick = 1:size(curDesign,2);
ax1.Box = 'off';
ax1.TickDir = 'out';
ax1.XTickLabel = {'sin' 'cos' 'constant' 'drift' 'x' 'y' 'z' 'pitch' 'roll' 'yaw' 'x''' 'y''' 'z''' 'pitch''' 'roll''' 'yaw'''};
ax1.XTickLabelRotation = -90;
ax1.YAxis.Color = 'none';
ax1.YAxis.Label.Visible = 'on';
ax1.YAxis.Label.Color = 'k';

%%
f = [f {figure('WindowStyle','docked')}];
curDesign = GLMs.random.hrDesign;
imagesc(curDesign,[-1 1]); colormap gray
title('Model-Free Response');
ax2 = gca; ax2.XTick = []; ax2.YTick = [];
ax2.XTick = 1:size(curDesign,2);
ax2.Box = 'off';
ax2.TickDir = 'out';
tLabel = cellstr([repmat('t_0+',11,1) num2str((1:11)','%-d')]);
ax2.XTickLabel = cat(1,tLabel,{'t_0base'},{'drift'});
ax2.XTickLabelRotation = -90;
ax2.YAxis.Color = 'none';
ax2.YAxis.Label.Visible = 'on';
ax2.YAxis.Label.Color = 'k';

% ax = [ax1 ax2];
% for i = 1:length(ax)
%     curAxe = ax(i);
%     curAxeB = axes();
%     linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
%     curAxeB.Color = 'none'; curAxeB.Box = 'on';
%     curAxeB.XTick = []; curAxeB.YTick = [];
% end


%% Fixed Effect model
f = [f {figure('WindowStyle','docked')}];
curDesign = GLMs.fixed.sinDesign;
imagesc(curDesign,[-1 1]); colormap gray
title({'Fixed-Effect Model' '\fontsize{8}Sinusoidal Response'});
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
ylabel('TRs');
ax1.XTick = [1.5 3.5 5.5];
ax1.Box = 'off';
ax1.TickDir = 'out';
ax1.XTickLabel = {'ori1' 'ori2' 'plaid'};
ax1.XTickLabelRotation = -90;
ax1.YAxis.Color = 'none';
ax1.YAxis.Label.Visible = 'on';
ax1.YAxis.Label.Color = 'k';
% ax1.XTick = [ax1.XTick mean([6.5 size(curDesign,2)+0.5])];
% ax1.XTickLabel = [ax1.XTickLabel; {'constant+drift for each run'}];
% curAxeB.TickDir = 'out';


ax = [ax1];
for i = 1:length(ax)
    curAxe = ax(i);
    curAxeB = axes();
    linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
    curAxeB.Color = 'none'; curAxeB.Box = 'on';
    curAxeB.XTick = []; curAxeB.YTick = [];
end
curAxeB.XTick = mean([6.5 size(curDesign,2)+0.5]);
curAxeB.XTickLabel = {'constant+drift for each run'};
curAxeB.TickDir = 'out';



