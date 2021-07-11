function ax = plotIm(ax,im,cLim)
cbWidth = 0.03;

alphaData = ones(size(im));
if ~exist('cLim','var') || isempty(cLim)
    cLim = [min(im(:)) max(im(:))];
end
cb = flipud(linspace(cLim(1),cLim(2),size(im,1))');
cb = repmat(cb,1,round(size(im,1).*cbWidth));
im = cat(2,cb,im);
alphaData = cat(2,zeros(size(cb)),alphaData);
hIm = imagesc(im,cLim);
hIm.YData = fliplr(hIm.YData);
hIm.AlphaData = alphaData;
ax.YDir = 'normal';
ax.PlotBoxAspectRatio = [1+cbWidth size(im,1)/size(im,2) 1];
xticks([]); yticks([]);
colormap(ax,'gray');
ax.Color = 'none';
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';