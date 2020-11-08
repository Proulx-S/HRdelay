function h = plotSVM_patAct7(svm,h,plotIt)

% if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
%     plotIt = svm.plotIt;
% else
%     plotIt = 1;
% end
if ~exist('plotIt','var')
    plotIt = 1;
end

%% Compute data to plot

patStr1_list = {'cart' 'pol'};
patStr2_list = {'pat1' 'pat2'};
depVar_list = {'amp' 'delay'};

hCount = 0;
for i = 1:length(patStr1_list)
    for ii = 1:length(patStr2_list)
        for iii = 1:length(depVar_list)
            hCount = hCount+1;
            
            
            patStr1 = patStr1_list{i};
            patStr2 = patStr2_list{ii};
            depVar = depVar_list{iii};
            curData = svm.r.crossInfo.f.(patStr1).(patStr2).(depVar);
            
            tmp1 = mean(curData(1:end/2,:),1);
            tmp2 = mean(curData(end/2+1:end,:),1);
            allMean.cart.pat1.amp = [tmp1' tmp2'];
            tmp1 = std(curData(1:end/2,:),[],1);
            tmp2 = std(curData(end/2+1:end,:),[],1);
            allStd.cart.pat1.amp = [tmp1' tmp2'];
            
            
            h(hCount) = figure('WindowStyle','docked');
            bar(allMean.cart.pat1.amp)
            set(gca,'xtickLabel',svm.r.crossInfo_dim1Train)
            xlabel('trained on')
            legend({'tested on grat1' 'tested on grat2'},'location','bestoutside')
            hold on
            xtick = get(gca,'xtick');
            offset = 0.14;
            errorbar(xtick-offset,allMean.cart.pat1.amp(:,1),allStd.cart.pat1.amp(:,1),'.')
            errorbar(xtick+offset,allMean.cart.pat1.amp(:,2),allStd.cart.pat1.amp(:,2),'.')
            title(patStr2)
            ylabel([patStr1 'pattern average of ' depVar ' (std)'])
        end
    end
end