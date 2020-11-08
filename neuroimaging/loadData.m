function data = loadData(dataRepo,dataDir,dataLevel,fileList)

switch dataLevel
    case 'z0'
        fileList
    case 'z'
        for sessInd = 1:size(fileList,2)
            for subjInd = 1:size(fileList,1)
                load(fullfile(dataRepo,dataDir,dataLevel,[fileList{subjInd,sessInd} '.mat']),'d')
                ori1{subjInd,sessInd} = d.xData(d.label==1,:);
                ori2{subjInd,sessInd} = d.xData(d.label==2,:);
                plaid{subjInd,sessInd} = d.normData;
            end
        end
        data.ori1 = ori1;
        data.ori2 = ori2;
        data.plaid = plaid;
end




