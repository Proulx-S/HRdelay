function polyCont = credibleInt2D(d,alpha)
warning('off','MATLAB:polyshape:repairedBySimplify')
nGridPts = 100;

xLim = [min(d(:,1)) max(d(:,1))];
yLim = [min(d(:,2)) max(d(:,2))];
deltaGrid = max([range(xLim) range(yLim)])/nGridPts;
x = (xLim(1)-0.15*range(xLim)):deltaGrid:(xLim(2)+0.15*range(xLim));
y = (yLim(1)-0.15*range(yLim)):deltaGrid:(yLim(2)+0.15*range(yLim));
[X,Y] = meshgrid(x,y);

[densityXY,~] = ksdensity(d,[X(:) Y(:)]);
density = nan(size(X));
density(:) = densityXY;
% imagesc(X(1,:),Y(:,1),density)
% contour(X,Y,density)


% Staircase down to the density threshold that gives the 2D interval
% closest to the specified alpha
level = sum(density(:))/numel(density)*[alpha alpha];
step = level*10;
delta = inf;
while abs(delta(end))>0.001
    % find resample points that are within the current above-threshold area
    M = contourc(x,y,density,level);
    polyCont = polyshape;
    while ~isempty(M)
        polyCont = addboundary(polyCont,M(1,2:1+M(2,1)),M(2,2:1+M(2,1)));
        M(:,1:1+M(2,1)) = [];
    end
    in = polyCont.isinterior(d(:,1),d(:,2));
    
    % compute the current p value
    curP = (1-nnz(in)/length(in));
    % compare to alpha level
    delta(end+1) = alpha - curP;
    
    % stoping rule
    [a,b] = unique(abs(delta));
    [~,c] = sort(abs(a));
    tmp = min(length(c),3);
    if sum(b(c(1:tmp)))/sum(b)>0.5
        break
    end
    
    % descend the gradient
    if diff(sign(delta(end-1:end)))
        step = step/2;
    end
    level = level + sign(delta(end))*step;
end

% imagesc(density); hold on
% contour(x,y,density,level);