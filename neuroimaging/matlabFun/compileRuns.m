function metricCompiled = compileRuns(metric,label,fitType,dataType)



%Prepare stuff and labels
selectedFields = fields(metric.(dataType){1}.(fitType));
% selectedFields = getFieldToCompile(metric);
nVox = size(metric.(dataType){1}.ind,1);
labelVal = unique(label)';
for i = 1:length(labelVal)
    labelName{i,1} = ['label' num2str(labelVal(i))];
    labelInd{i,1} = find(labelVal(i)==label);
end
tmp = label; clear label;
label.label = tmp;
label.labelName = labelName;
label.labelVal = labelVal;
label.labelInd = labelInd;
clear labelName labelVal labelInd
nLabel = length(label.labelName);
nRun = length(metric.(dataType));
nRunPerLabel = nRun/nLabel;
obsDim = find(nVox~=size(metric.(dataType){1}.(fitType).(selectedFields{1})));
nObsPerRun = size(metric.(dataType){1}.(fitType).(selectedFields{1}),obsDim);


nObsPerLabel = length(metric.(dataType))*nObsPerRun/nLabel;



%Compile data
for curFieldInd = 1:length(selectedFields)
    if strcmp(selectedFields(curFieldInd),'RMSEall')
        dim4Size = size(metric.(dataType){1}.(fitType).(selectedFields{curFieldInd}),2);
        dim5Size = size(metric.(dataType){1}.(fitType).(selectedFields{curFieldInd}),3);
        curDataAll = nan(nVox,nObsPerLabel,nLabel,dim4Size,dim5Size);
    else
        curDataAll = nan(nVox,nObsPerLabel,nLabel);
    end
    curSessionAll = nan(1,nObsPerLabel,nLabel);
%     curData = nan(nVox,nObsPerRun,1);
%     curDataSession = nan(nVox,nObsPerLabel,nLabel);
    for curRun = 1:nRun
%         if strcmp(selectedFields(curFieldInd),'RMSEall')
%         else
%             curData = metric.(dataType){curRun}.(fitType).(selectedFields{curFieldInd});
%             switch obsDim
%                 case 1
%                     curData = metric.(dataType){curRun}.(fitType).(selectedFields{curFieldInd})';
%                 case 2
                    curData = metric.(dataType){curRun}.(fitType).(selectedFields{curFieldInd});
%             end
%         end
        
        curSession = num2str(metric.runList(curRun));
        curSession = repmat(str2double(curSession(1)),1,nObsPerRun);
        
        curLabelInd = find(label.label(curRun)==label.labelVal);
        curRunInd = find(label.labelInd{curLabelInd}==curRun);
        curObsInd = ((curRunInd-1)*nObsPerRun+1):(curRunInd*nObsPerRun);
        if strcmp(selectedFields(curFieldInd),'RMSEall')
            curDataAll(:,curObsInd,curLabelInd,:,:) = reshape(curData,[size(curData,1) 1 1 size(curData,2) size(curData,3)]);
        else
            curDataAll(:,curObsInd,curLabelInd) = curData;
        end
        curSessionAll(1,curObsInd,curLabelInd) = curSession;
    end
        
%         
%         for obs = 1:nObsPerLabel
%             
%             curRunInd = label.labelInd{curLabelInd}(obs);
%             curData(:,obs,curLabelInd) = metric.(dataType){curRun}.(fitType).(selectedFields{curFieldInd});
%         end
%     end
%         
%         
%         
%     for obs = 1:nObsPerLabel
%         for curLabelInd = 1:nLabel
%             curRunInd= label.labelInd{curLabelInd}(obs);
%             curData(:,obs,curLabelInd) = metric.(dataType){curRunInd}.(fitType).(selectedFields{curFieldInd});
%         end
%     end
    metricCompiled.(selectedFields{curFieldInd}) = curDataAll;
end
metricCompiled.label = label;
metricCompiled.voxInd = metric.(dataType){1}.ind;
metricCompiled.fittingParam = metric.(dataType){1}.p;
metricCompiled.dataFields = selectedFields;
metricCompiled.sessionLabel = curSessionAll;


% %Get session labels
% curDataSession = nan(1,nObsPerLabel,nLabel);
% for obs = 1:nObsPerLabel
%     for curLabelInd = 1:nLabel
%         curRunInd = label.labelInd{curLabelInd}(obs);
%         tmp = num2str(metric.runList(curRunInd));
%         curDataSession(1,obs,curLabelInd) = str2num(tmp(1)); clear tmp
%     end
% end
% metricCompiled.sessionLabel = curDataSession;







% function selectedFields = getFieldToCompile(metric)
% fields = fieldnames(metric.(dataType){1});
% selectedFields = {};
% i = 0;
% done = false;
% while ~done
%     i = i+1;
%     curField = fields{i};
%     curFieldData = metric.(dataType){1}.(curField);
%     if isnumeric(curFieldData);
%         if i==1
%             targetSize = size(curFieldData);
%         else
%             if size(curFieldData)~=targetSize
%                 done = true;
%             end
%         end
%     else
%         done = true;
%     end
%     if ~done
%         selectedFields = [selectedFields; {curField}];
%     end
% end
