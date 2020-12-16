function plotDecodingPerm(res,saveFig)
if ~exist('saveFig','var') || isempty(saveFig)
    saveFig = 0;
end

spaceList = fields(res)';
for spaceInd = 1:length(spaceList)
    
    res.(spaceList{spaceInd}).perm.acc = permute(res.(spaceList{spaceInd}).perm.acc,[2 3 1]);
    
    % Group
    fGroup = figure('WindowStyle','docked');
    acc = res.(spaceList{spaceInd}).perm.summary.acc;
    nObs = sum(res.(spaceList{spaceInd}).nObs(:));
    hit = round(acc.*nObs);
    plotHitsAndBino(hit,nObs);
    xlabel('hit count')
    title(spaceList{spaceInd},'interpreter','none')
    
    if saveFig
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,[spaceList{spaceInd} '__group']);
        fGroup.Color = 'none';
        set(findobj(fGroup.Children,'type','Axes'),'color','none')
        saveas(fGroup,[filename '.svg']); disp([filename '.svg'])
        fGroup.Color = 'w';
        saveas(fGroup,filename); disp([filename '.fig'])
        saveas(fGroup,filename); disp([filename '.jpg'])
    end


    % Individual subjects
    accAll = res.(spaceList{spaceInd}).perm.acc;
    nObsAll = res.(spaceList{spaceInd}).nObs;
    hitAll = round(accAll.*nObsAll);
    yLim = [];
    fSubj = figure('WindowStyle','docked');
    for subjInd = 1:size(accAll,1)
        for sessInd = 1:size(accAll,2)
            subplot(size(accAll,1),size(accAll,2),sessInd + (subjInd-1)*2)
            nObs = nObsAll(subjInd,sessInd);
            hit = squeeze(hitAll(subjInd,sessInd,:));
            plotHitsAndBino(hit,nObs);
            xticks([0 nObs/2 nObs]);
            xlabel(''); ylabel('')
            yLim = [yLim; ylim];
        end
    end
    yLim = [min(yLim(:,1)) max(yLim(:,2))];
    for subjInd = 1:size(accAll,1)
        for sessInd = 1:size(accAll,2)
            subplot(size(accAll,1),size(accAll,2),sessInd + (subjInd-1)*2)
            ylim(yLim);
            switch subjInd
                case 1
                    title(['sess' num2str(sessInd)])
                case size(accAll,1)
                    xlabel('Hit count')
            end
            if sessInd==1
                ylabel([res.(spaceList{spaceInd}).subjList(subjInd) {'Density'}])
            end
        end
    end
    
    if saveFig
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,[spaceList{spaceInd} '__subj']);
        fGroup.Color = 'none';
        set(findobj(fGroup.Children,'type','Axes'),'color','none')
        saveas(fGroup,[filename '.svg']); disp([filename '.svg'])
        fGroup.Color = 'w';
        saveas(fGroup,filename); disp([filename '.fig'])
        saveas(fGroup,filename); disp([filename '.jpg'])
    end
end


function pd1 = plotHitsAndBino(hit,n)
% --- Plot data originally in dataset "hit data"
[CdfF,CdfX] = ecdf(hit,'Function','cdf');  % compute empirical cdf
BinEdge = (0:n+1)-0.5;
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist'); hold on
hLine.FaceColor = 'k';
xlim([0 n])
xlabel('Hit count');
ylabel('Density')

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XIncr = max(1,floor((XLim(2)-XLim(1))/100));
XGrid = floor(XLim(1)):XIncr:ceil(XLim(2));


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('binomial',[ 142, 0.5])
pd1 = fitdist(hit, 'binomial', 'n', n);
YPlot = pdf(pd1,XGrid);
plot(XGrid,YPlot,'Color','r');
