function data = avVox(data)
for i = 1:numel(data.ori1)
    data.ori1{i} = mean(data.ori1{i},2);
    data.ori2{i} = mean(data.ori2{i},2);
    data.plaid{i} = mean(data.plaid{i},2);
end
% for sessInd = 1:size(fileList,2)
%     for subjInd = 1:size(fileList,1)
%         data.ori1{subjInd,sessInd} = mean(data.ori1{subjInd,sessInd},2);
%         data.ori2{subjInd,sessInd} = mean(data.ori2{subjInd,sessInd},2);
%         data.plaid{subjInd,sessInd} = mean(data.plaid{subjInd,sessInd},2);
%     end
% end
