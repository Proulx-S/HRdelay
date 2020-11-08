function metricOut = applyThresh(metricIn,treshMetric,threshField,thresh,method)
if ~exist('method','var')
    method = 'acrossRunAverage_anyGrat';
end

metricOut = metricIn;

threshData = treshMetric.(threshField);
if strcmp(threshField,'p')
    threshData = 1-threshData;
end

switch method
    case 'acrossRunAverage_anyGrat'
        threshLabelsInd = [1 2];
        threshData = threshData(:,:,threshLabelsInd);

        threshData = mean(threshData,2);
        aboveThreshInd = threshData>=thresh;
        aboveThreshInd = any(aboveThreshInd,3);
    case 'acrossRunAverage_anyGratOrPlaid'
        threshLabelsInd = [1 2 3];
        threshData = threshData(:,:,threshLabelsInd);

        threshData = mean(threshData,2);
        aboveThreshInd = threshData>=thresh;
        aboveThreshInd = any(aboveThreshInd,3);
end





for curFieldInd = 1:length(metricIn.dataFields)
    metricOut.(metricIn.dataFields{curFieldInd}) = metricIn.(metricIn.dataFields{curFieldInd})(aboveThreshInd,:,:);
end
metricOut.voxInd = metricIn.voxInd(aboveThreshInd,:);
