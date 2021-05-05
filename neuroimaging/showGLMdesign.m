function showGLMdesign(GLMs)

f = figure('WindowStyle','docked');

subplot(1,16+25,25+(1:16))
curDesign = GLMs.random.sinDesign;
imagesc(curDesign,[-1 1]); colormap gray
title('Sinusoidal Response');
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
%     xlabel('Regressors');
ylabel('TRs');
ax1.XTick = 1:size(curDesign,2);
ax1.Box = 'off';
ax1.TickDir = 'out';
ax1.XTickLabel = {'sin' 'cos' 'constant' 'drift' 'x' 'y' 'z' 'pitch' 'roll' 'yaw' 'x''' 'y''' 'z''' 'pitch''' 'roll''' 'yaw'''};
ax1.XTickLabelRotation = -90;
ax1.YAxis.Color = 'none';
ax1.YAxis.Label.Visible = 'on';
ax1.YAxis.Label.Color = 'k';

xVal = diff(ax1.XLim);
xWidth = ax1.Position(3);

subplot(1,16+25,0+(1:25))
% GLMs.hrDesign(:,1:end-2)
curDesign = GLMs.random.hrDesign;
imagesc(curDesign,[-1 1]); colormap gray
title('Model-Free Response');
ax2 = gca; ax2.XTick = []; ax2.YTick = [];
%     xlabel('Regressors');% ylabel('TRs');
ax2.XTick = 1:size(curDesign,2);
ax2.Box = 'off';
ax2.TickDir = 'out';
tLabel = cellstr([repmat('t_0+',11,1) num2str((1:11)','%-d')]);
ax2.XTickLabel = cat(1,tLabel,{'t0 (base)'},{'drift'});
ax2.XTickLabelRotation = -90;
ax2.YAxis.Color = 'none';
ax2.YAxis.Label.Visible = 'on';
ax2.YAxis.Label.Color = 'k';

suptitle('Design Matrices (random-effect)');

ax1.XAxis.FontSize = ax1.XAxis.FontSize*0.8;
ax2.XAxis.FontSize = ax2.XAxis.FontSize*0.8;

ax = [ax1 ax2];
for i = 1:length(ax)
    curAxe = ax(i);
    curAxeB = axes();
    linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
    curAxeB.Color = 'none'; curAxeB.Box = 'on';
    curAxeB.XTick = []; curAxeB.YTick = [];
end


% Fixed Effect model
f = figure('WindowStyle','docked');
curDesign = GLMs.fixed.sinDesign;
imagesc(curDesign,[-1 1]); colormap gray


run1 = results.OLS.fixed.designmatrix(sessInd,:);
runOrderLabel = results.OLS.fixed.designmatrixRunInd(sessInd,:);
run1 = run1(:,~all(run1==0,1));
runOrderLabelList = unique(runOrderLabel);
tmp = [];
for i = 1:length(runOrderLabelList)
    trInd = runOrderLabel==runOrderLabelList(i);
    tmp = cat(1,tmp,run1(trInd,:));
end
tmp(:,7:end) = run1(:,7:end);
run1 = tmp;
tmp = run1(:,7:2:end);
run1(:,7:2:end) = run1(:,7:2:end)./mode(tmp(tmp~=0));
tmp = run1(:,8:2:end);
run1(:,8:2:end) = run1(:,8:2:end)./max(abs(tmp(tmp~=0)));

imagesc(run1); colormap gray
title({'Fixed-Effect Model' '\fontsize{8}Sinusoidal Response'});
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
%     xlabel('Regressors');
ylabel('TRs');
ax1.XTick = [1.5 3.5 5.5];
ax1.Box = 'off';
ax1.TickDir = 'out';
ax1.XTickLabel = {'ori1' 'ori2' 'plaid'};
ax1.XTickLabelRotation = -90;
ax1.YAxis.Color = 'none';
ax1.YAxis.Label.Visible = 'on';
ax1.YAxis.Label.Color = 'k';

ax = [ax1];
for i = 1:length(ax)
    curAxe = ax(i);
    curAxeB = axes();
    linkprop([curAxe curAxeB],'position'); linkprop([curAxe curAxeB],'XLim');
    curAxeB.Color = 'none'; curAxeB.Box = 'on';
    curAxeB.XTick = []; curAxeB.YTick = [];
    %         delta = mode(diff(curAxe.XTick));
    %         curAxeB.XTick = [curAxe.XTick curAxe.XTick(end)+delta]-delta/2;
    %         curAxeB.TickLength = [0 0];
    %         curAxeB.XTickLabel = [];
    %         curAxeB.XGrid = 'on';
    %         curAxeB.GridAlpha = 1;
end
curAxeB.XTick = mean([6.5 size(run1,2)+0.5]);
curAxeB.XTickLabel = {'constant+drift for each run'};
curAxeB.TickDir = 'out';

if saveFig
    saveas(f,fullfile(repo,funDir,outDir,'designMatrices_fixedEffect'))
    if verbose
        disp(fullfile(repo,funDir,outDir,'designMatrices_fixedEffect.fig'))
    end
end


