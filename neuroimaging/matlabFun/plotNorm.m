function f0 = plotNorm(d,p,featSel,voxFlag,visibilityFlag)
if ~exist('voxFlag','var') || isempty(voxFlag)
    voxFlag = 0;
end
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

%% Initate stuff
featSelConds_labelList = featSel.featSeq.condPairList;
for i = 1:length(featSelConds_labelList)
    featSelConds_labelList{i} = num2str(featSelConds_labelList{i});
end
ind_nSpecFeatSelCond = [1 2 3];
ind_nSpecFeatSelCond = squeeze(ismember(featSelConds_labelList,num2str(ind_nSpecFeatSelCond)));

featSelList = featSel.featSeq.featSelList;
for featInd = 1:length(featSelList)
    tmp = strsplit(featSel.featSeq.featSelList{featInd},': ');
    featSelList(featInd) = tmp(1);
end

%% Get data
[X,y,~] = getXYK(d,p);

%% Plot
% select most active voxel
% featInd = ismember(featSelList,'act');
featInd = ismember(featSelList,'respVecSig');
% featInd = ismember(featSelList,'respVecDiff');
[~,b] = sort(featSel.featSeq.featVal(featSel.indIn,featInd,ind_nSpecFeatSelCond),'descend');
f0 = plotPolNormExample(X(:,featSel.indIn),y,p,b(1),voxFlag,visibilityFlag);
% f1 = plotPolNormExampleVox(X,SVMspace,b(1));
% f2 = plotPolNormExampleRep(X,y,SVMspace,b(1));

% x = cat(3,X(y==1,b(1)),X(y==2,b(1)));
% f2 = plotPolNormExample(x,SVMspace);


function f = plotPolNormExample(x,y,p,b,voxFlag,visibilityFlag)
if ~exist('voxFlag','var')
    voxFlag = 0;
end
if ~exist('visibilityFlag','var')
    visibilityFlag = 1;
end

if visibilityFlag
    f = figure('WindowStyle','docked','visible','on');
else
    f = figure('WindowStyle','docked','visible','off');
end

%% Polar Normalization
% Plot before
subplot(2,2,1); clear hPP
if voxFlag
    hPP3 = polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
    hLeg = legend(hPP3,{'voxels'},'box','on');
else
    hPP1 = polarplot(angle(x(y==1,b)),abs(x(y==1,b)),'.'); hold on
    hPP2 = polarplot(angle(x(y==2,b)),abs(x(y==2,b)),'.'); hold on
    hPP3 = polarplot(angle(x(1,:)),abs(x(1,:)),'.k'); hold on
    
    uistack(hPP3,'bottom');
    hPP1.MarkerSize = hPP1.MarkerSize*2;
    hPP2.MarkerSize = hPP2.MarkerSize*2;
    hLeg = legend([hPP1 hPP2 hPP3],{'1vox; Areps' '1vox; Breps' 'allVox; 1rep'},'box','on');
end
hPP3.MarkerSize = eps;
hPP3.Color = [1 1 1].*0;
hLeg.Location = 'southeast';
drawnow
title('before polarNorm')
drawnow

if voxFlag
    subplot(2,2,3);
    histogram(abs(mean(x,1)))
    xlabel('pre-norm amp (%BOLD)')
    subplot(2,2,4);
    meanDelay = angle(mean(x(:)));
end

% Normalize
x = polarSpaceNormalization(x,p.chanSpace);

if voxFlag
    subplot(2,2,4);
    tmp = -(angle(mean(x,1)) + meanDelay)./pi.*6;
    histogram(tmp)
    xlabel('pre-norm delay (sec)')
end


% Plot after
subplot(2,2,2);
if voxFlag
    hPP3 = polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
else
    hPP1 = polarplot(angle(x(y==1,b)),abs(x(y==1,b)),'.'); hold on
    hPP2 = polarplot(angle(x(y==2,b)),abs(x(y==2,b)),'.'); hold on
    hPP3 = polarplot(angle(x(1,:)),abs(x(1,:)),'.k'); hold on
    uistack(hPP3,'bottom');
    hPP1.MarkerSize = hPP1.MarkerSize*2;
    hPP2.MarkerSize = hPP2.MarkerSize*2;
end
hPP3.MarkerSize = eps;
hPP3.Color = [1 1 1].*0;
title('after polarNorm')
ax = gca;
ax.RLim = [0 prctile(hPP3.RData,95)];
drawnow


if ~voxFlag
    %% Cartesian Normalization
    % Plot before
    subplot(2,2,3);
    hScat1 = scatter(real(x(y==1,b)),imag(x(y==1,b)),'o','filled'); hold on
    hScat2 = scatter(real(x(y==2,b)),imag(x(y==2,b)),'o','filled'); hold on
    hScat3 = scatter(real(x(1,:)),imag(x(1,:)),'ko','filled'); hold on
    % hScat3 = scatter(real(x(:)),imag(x(:)),'ko','filled'); hold on
    % hScat3 = scatter(real(x(1,:)),imag(x(1,:)),'ko','filled'); hold on
    uistack(hScat3,'bottom')
    hScat3.SizeData = hScat3.SizeData./8;
    hScat1.MarkerEdgeColor = 'w';
    hScat2.MarkerEdgeColor = 'w';
    ax = gca;
    ax.PlotBoxAspectRatio = [1 1 1];
    % lim = prctile(abs([real(x(1,:)) imag(x(1,:))]),95);
    xLim = prctile(real(x(1,:)),[2.5 97.5]);
    yLim = prctile(imag(x(1,:)),[2.5 97.5]);
    % lim = [-lim lim];
    if ~all(xLim==[0 0])
        xlim(xLim);
    end
    if ~all(yLim==[0 0])
        ylim(yLim);
    end
    % switch p.chanSpace
    %     case {'cart' 'cartReal'}
    %         xlim(lim);
    %         ylim(lim);
    %     case 'cartNoAmp'
    %     case 'cartNoDelay'
    %         xlim([0 lim(2)]);
    %     otherwise
    %         error('X')
    % end
    % xLim = xlim;
    delta = abs(diff(xLim)).*0.1;
    if ~(xLim(1)<0)
        xLim(1) = -delta;
    end
    if ~(xLim(2)>0)
        xLim(2) = +delta;
    end
    xlim(xLim)
    
    yLim = ylim;
    if ~(yLim(1)<0)
        yLim(1) = -delta;
    end
    if ~(yLim(2)>0)
        yLim(2) = +delta;
    end
    ylim(yLim)
    
    uistack(plot([0 0],ylim,'-k'),'bottom');
    uistack(plot(xlim,[0 0],'-k'),'bottom');
    grid on
    title('before cartNorm')
    xlabel('real')
    ylabel('imag')
    drawnow
    hLeg = legend([hScat1 hScat2 hScat3],{'1vox; Breps' '1vox; Breps' 'allVox; 1rep'},'box','on');
    hLeg.Location = 'northwest';
    drawnow
    
    % Normalize
    x = cartSpaceNormalization(x,p.chanSpace);
    
    %Plot after
    subplot(2,2,4);
    hScat1 = scatter(real(x(y==1,b)),imag(x(y==1,b)),'o','filled'); hold on
    hScat2 = scatter(real(x(y==2,b)),imag(x(y==2,b)),'o','filled'); hold on
    hScat3 = scatter(real(x(1,:)),imag(x(1,:)),'ko','filled'); hold on
    % hScat3 = scatter(real(x(:)),imag(x(:)),'ko','filled'); hold on
    % hScat3 = scatter(real(x(1,:)),imag(x(1,:)),'ko','filled'); hold on
    uistack(hScat3,'bottom')
    hScat3.SizeData = hScat3.SizeData./8;
    hScat1.MarkerEdgeColor = 'w';
    hScat2.MarkerEdgeColor = 'w';
    ax = gca;
    ax.DataAspectRatio = [1 1 1];
    ax.PlotBoxAspectRatio = [1 1 1];
    xLim = prctile(real(x(1,:)),[2.5 97.5]);
    yLim = prctile(imag(x(1,:)),[2.5 97.5]);
    % lim = prctile(abs([real(x(:)) imag(x(:))]),95);
    % lim = [-lim lim];
    if ~all(xLim==[0 0])
        xlim(xLim);
    end
    if ~all(yLim==[0 0])
        ylim(yLim);
    end
    % switch p.chanSpace
    %     case {'cart' 'cartReal'}
    %     case 'cartNoAmp'
    %     case 'cartNoDelay'
    %         xlim([0 lim(2)]);
    %     otherwise
    %         error('X')
    % end
    xLim = xlim;
    delta = abs(diff(xLim)).*0.1;
    if ~(xLim(1)<0)
        xLim(1) = -delta;
    end
    if ~(xLim(2)>0)
        xLim(2) = +delta;
    end
    xlim(xLim)
    
    yLim = ylim;
    delta = abs(diff(yLim)).*0.1;
    if ~(yLim(1)<0)
        yLim(1) = -delta;
    end
    if ~(yLim(2)>0)
        yLim(2) = +delta;
    end
    ylim(yLim)
    
    uistack(plot([0 0],ylim,'-k'),'bottom');
    uistack(plot(xlim,[0 0],'-k'),'bottom');
    grid on
    title('after cartNorm')
    xlabel('real')
    ylabel('imag')
    drawnow
end

