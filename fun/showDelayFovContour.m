function f = showDelayFovContour(X,Y,vecXY,cMap,U,V,cont,contInd,visibilityFlag)
if ~exist('visibilityFlag','var') || isempty(visibilityFlag)
    visibilityFlag = 'on';
end
if ~exist('contInd','var') || isempty(contInd)
    contInd = 1:length(cont);
end

f = figure('WindowStyle','docked','visible',visibilityFlag);
% [ha, ~] = tight_subplot(2,1,0,0.03,0.03); drawnow
%% Get contour on smoothed map
% set(f,'CurrentAxes',ha(1))
if ~isempty(X)
    imagesc(X(1,:),Y(:,1),vecXY,[0 pi]); hold on
end
ax = gca;
colormap(flip(cMap,1));
set(ax,'YDir','normal')
hScat = scatter(U,V); hold on
hScat.MarkerEdgeColor = 'none';
hScat.MarkerFaceColor = 'k';
hScat.SizeData = 4^2;
hScat.MarkerFaceAlpha = 0.1;
ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
% [M,c] = contour(X,Y,vecXYcont,ones(2,1).*level); hold on

for i = 1:length(cont)
    if ismember(i,abs(contInd)) | contInd==inf
        plot(cont{i}(:,1),cont{i}(:,2),'k')
    end
end

