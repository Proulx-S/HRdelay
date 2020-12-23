function accPerm = plotDecodingPerm_acc(res,figOption)
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
figOption.subj = 1; % 'all' or subjInd

spaceList = fields(res)';
for spaceInd = 1:length(spaceList)
    nPerm = size(res.(spaceList{spaceInd}).perm.acc,3);
    
    % Group
    fGroup = figure('WindowStyle','docked');
    acc = res.(spaceList{spaceInd}).perm.summary.acc;
    nObs = sum(res.(spaceList{spaceInd}).nObs(:));
    hit = round(acc.*nObs);
    plotHitsAndBino(hit,nObs);
    xlabel('Hit rate')
    xlabel({'Hit rate' '\fontsize{8}(all subj and sess concatenated)'})
    ylabel({'Density' ['\fontsize{8}' num2str(nPerm) ' permutations']})
    title(spaceList{spaceInd},'interpreter','none')
    ax = gca;
    xTick = [0 nObs.*0.5 nObs];
    ax.XTick = xTick;
    ax.XTickLabel = compose(['%d/' num2str(nObs)],xTick);
    ax.XAxis.TickDirection = 'out';
    ax.Box = 'off';
    
    if figOption.save
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
    if figOption.subj
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
                ax = gca;
                xTick = [0 nObs.*0.5 nObs];
                ax.XTick = xTick;
                ax.XTickLabel = compose(['%d/' num2str(nObs)],xTick);
                ax.XAxis.TickDirection = 'out';
                ax.Box = 'off';
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
                        xlabel('Hit rate')
                        if sessInd==1
                            ylabel({'Density' ['\fontsize{8}' num2str(nPerm) ' permutations']})
                        end
                end
                if sessInd==2
                    h = ylabel(res.(spaceList{spaceInd}).subjList(subjInd))
                    ax = gca;
                    ax.YAxisLocation = 'right';
                    ax.YTick = [];
                    ax.YAxis.Visible = 'off';
                    ax.YAxis.Label.Visible = 'on';
                    ax.YAxis.Label.Rotation = 0;
                    ax.YAxis.Label.HorizontalAlignment = 'left';
                    ax.YAxis.Label.VerticalAlignment = 'middle';
                end
            end
        end
        h = suptitle(spaceList{spaceInd});
        h.Interpreter = 'none';
        
        if figOption.save
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
    accPerm.(spaceList{spaceInd}) = acc;
end



