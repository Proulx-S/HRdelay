function out = runSVM_patAct2(d,p,verbose,rep,doPerm)

if ~exist('verbose','var')
    verbose = 0;
end
if ~isfield(p,'rotTo')
    if strcmp(p.keepInfo,'i')
        p.rotTo = 0;
    else
        p.rotTo = pi/4;
    end
end
if ~exist('doPerm','var')
    doPerm = 0;
end

%% Permute data if needed
if doPerm
    d.xData = d.xData(randperm(size(d.xData,1)),:);
end

%% Define k-folding
if ischar(p.k)
    switch p.k
        %         case {'auto','autoCmplt'}
        %             k = length(d.crossVal)/2;
        case 'autoRun'
            d.crossVal = defineCrossVal(d,p);
            k = length(d.crossVal)/2/str2double(p.split(1));
        otherwise
            error('');
    end
else
    error('')
    k = p.k;
    p.k = 'randK';
end


%% Loop over cross-validation folds
dOrig = d;
out.hitRate = nan(k,1);
out.hitRate_d = nan(k,1);
out.hitRate_pat = nan(k,1);
out.hitRate_patRel = nan(k,1);

out.pat = nan(p.nObs,2);
out.patRel = nan(p.nObs,2);
keyboard
for fold = 1:k
    
    [tr,te,teInd1,teInd2] = prepareDataAtFold(dOrig,p,'r',fold,k,rep,verbose);
    
%     %% Sorting for feature slection (previously done with newFeatSorting.m)
%     % compute sorting
%     if ~exist('sortIndFix','var')
%         switch p.keepInfo
%             case {'m','r','i'}
%                 p.sort.feat = computeMI(tr.amp,tr.label==1,tr.label==2)';
%             case 'p'
%                 for dim = 1:size(tr.delay,2)
%                     [pval(dim), table] = circ_wwtest(tr.delay(tr.label==1,dim),tr.delay(tr.label==2,dim),[],[],0);
%                 end
%                 p.sort.feat = 1-pval;
%             case 'mp'
%                 [X,Y] = pol2cart(tr.delay,tr.amp);
%                 p.sort.feat = computeMI(complex(X,Y),tr.label==1,tr.label==2)';
%             otherwise
%                 error('did not code that')
%         end
%     else
%         p.sort.feat = nan(1,p.nFeatures);
%     end
    
%     % Sorting for functional ROI
%     p.sort.fun = p.funcROI.vec.stats;
%     [~,p.sort.indFun] = sort(p.sort.fun,'descend');
%     
%     trPreFunSort = tr;
%     tePreFunSort = te;
    
    keyboard
    % Run SVM
    %Train
    [out] = finalPrep(tr,p,'r');
    
    
    %Train and test
    outk = processDataAtFold(tr,te,p,'r');
    allFields = fields(out);
    for i = 1:length(allFields)
        if strfind(allFields{i},'hitRate')
            out.(allFields{i})(fold) = outk.(allFields{i});
        elseif strfind(allFields{i},'pat')
            out.(allFields{i})([teInd1;teInd2],:) = outk.(allFields{i});
        else
            error('something wrong here')
        end
    end
    clear outk
    
%     %% Loop over functional ROI selection levels
%     for funLevInd = 1:length(p.featLevList)
% %         % Apply fun sorting (and select fun)
% %         trFeatSort.amp = trPreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
% %         trFeatSort.delay = trPreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
% %         teFeatSort.amp = tePreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
% %         teFeatSort.delay = tePreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
% %         curFun = p.sort.fun(p.sort.indFun((1:p.featLevList(funLevInd))));
% %         curFeat = p.sort.feat(p.sort.indFun((1:p.featLevList(funLevInd))));
% %         
% %         % Sort feat
% %         [~,curIndFeat] = sort(curFeat,'descend');
% %         
% %         % Apply feat sorting
% %         trFeatSort.amp = trFeatSort.amp(:,curIndFeat);
% %         trFeatSort.delay = trFeatSort.delay(:,curIndFeat);
% %         teFeatSort.amp = teFeatSort.amp(:,curIndFeat);
% %         teFeatSort.delay = teFeatSort.delay(:,curIndFeat);
% %         curFun = curFun(curIndFeat);
% %         curFeat = curFeat(curIndFeat);
%         
% 
% 
% 
% 
%             
%             
%         %% Loop over feature selection levels
%         for featLevInd = 1:length(p.featLevList)
%             if p.featLevList(featLevInd)>length(curFun)
%                 break
%             end
% %             % Select feat
% %             tr.amp = trFeatSort.amp(:,1:p.featLevList(featLevInd));
% %             tr.delay = trFeatSort.delay(:,1:p.featLevList(featLevInd));
% %             te.amp = teFeatSort.amp(:,1:p.featLevList(featLevInd));
% %             te.delay = teFeatSort.delay(:,1:p.featLevList(featLevInd));
% %             curFun(1:p.featLevList(featLevInd));
% %             curFeat(1:p.featLevList(featLevInd));
%             
% %             % Run SVM
% %             %Train and test
% %             out = processDataAtFold(tr,te,p,'r',fold);
%             
%             
%             
%         end
%     end
    %     if p.doSVM
    %         %Sort back weights to original and store
    %         w(:,:,fold) = sort_back(wFold,repmat(sortInd(fold,:),size(wFold,1),1),2);
    %         A(:,:,fold) = sort_back(Afold,repmat(sortInd(fold,:),size(Afold,1),1),2);
    %         %         data(:,:,fold) = sort_back(dataFold,repmat(sortInd(fold,:),size(dataFold,1),1),2);
    %     end
    if verbose
        toc
    end
    
    
end
% %Average over cross-validation folds
% switch p.algorithm
%     case 'runSVM'
%         out.hitRate = mean(out.hitRate,1);
%     case 'runSVM_RFE'
%         %         out.hitRate = mean(mean(mean(out.hitRate,1),2),4);
%         out.hitRate = mean(mean(out.hitRate,1),2);
%     case 'runSVM_patAct'
%         %         out.hitRate = mean(out.hitRate,1);
%         %         out.hitRate_d = mean(out.hitRate_d,1);
%         %         out.hitRate_pat = mean(out.hitRate_pat,1);
%         %         out.hitRate_patRel = mean(out.hitRate_patRel,1);
% end


