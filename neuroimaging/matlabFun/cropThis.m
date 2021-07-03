function [brain,lim] = cropThis(brain,lim)

if ~exist('lim','var') || isempty(lim)
    for label = {'X' 'Y'}
        figure('WindowStyle','docked')
        switch char(label)
            case 'X'
                brainX = squeeze(mean(brain,2));
            case 'Y'
                brainX = squeeze(mean(brain,1));
            otherwise
                error('X')
        end
        xLim = [1 length(brainX)];
        x = xLim(1):xLim(2);
        plot(x,brainX); hold on
        xlim(xLim)
        yLim = ylim;
        y = yLim(1):yLim(2);
        tmpIm = ones(diff(xLim)+1,diff(yLim)+1);
        tmpIm(1,1) = 0;
        hIm = imagesc(x,y,tmpIm); colormap gray
        uistack(hIm,'bottom')
        drawnow
        roipoly
        L1 = input('PolyCordinate = ');
        plot(L1(:,1),L1(:,2))
        x1 = x;
        xLim = nan(size(brainX,2),2);
        for i = 1:size(brainX,2)
            y1 = brainX(:,i)';
            x2 = L1(:,1)';
            y2 = L1(:,2)';
            P = InterX([x1;y1],[x2;y2]);
            plot(x2,y2,P(1,:),P(2,:),'ro')
            xLim(i,:) = P(1,:);
        end
        lim.(char(label)) = xLim;
    end
    plotFlag = 1;
else
    plotFlag = 0;
end


cropMask = false(size(brain));
if plotFlag
    figure('WindowStyle','docked')
    hTile = tiledlayout(size(brain,3),1,'TileIndexing','columnmajor');
end

for i = 1:size(brain,3)
    if plotFlag
        nexttile
    end
    x = floor(lim.X(i,1)):ceil(lim.X(i,2));
    y = floor(lim.Y(i,1)):ceil(lim.Y(i,2));
    cropMask(x,y,i) = true;
    tmp = brain(:,:,i);
    m = max(brain(:));
    tmp(~cropMask(:,:,i)) = m+eps;
    if plotFlag
        imagesc(tmp)
    end
    brain(:,:,i) = tmp;
    if plotFlag
        xticks([]); yticks([]);
        ax = gca;
        ax.PlotBoxAspectRatio = [1 1 1];
    end
end
if plotFlag
    colormap gray
    hTile.TileSpacing = 'tight';
end
