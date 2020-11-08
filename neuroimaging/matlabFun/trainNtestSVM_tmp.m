function [hitRate,svmStruct] = trainNtestSVM_tmp(tr,trLabel,te,teLabel,p)
global err
autoscaleFlag = true;

%% Train the model
switch p.kernel
    case 'polynomial'
        svmParams = statset('MaxIter',1500000,'Display','off');
%         try
            svmStruct = svmtrain(tr,trLabel,'kernel_function',p.kernel,'polyorder',p.polyorder,'autoscale',autoscaleFlag,'method','SMO','boxconstraint',p.C,'options',svmParams);
%         catch err
%             return
%         end
    otherwise
        svmParams = statset('MaxIter',1500000,'Display','off');
%         try
            svmStruct = svmtrain(tr,trLabel,'kernel_function',p.kernel,'autoscale',autoscaleFlag,'method','SMO','boxconstraint',p.C,'options',svmParams);
%         catch err
%             warning(err.identifier,err.message)
%             return
%         end
end

%% Test the model
predictedLabel = svmclassify(svmStruct,te);
hitRate = length(find(predictedLabel==teLabel))/length(teLabel);
