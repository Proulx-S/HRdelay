clear all
close all

main = load('C:\Users\sebas\Desktop\main.mat');
branch = load('C:\Users\sebas\Desktop\branch.mat');

ans = [];
for subjInd = 1:length(main.dAll)
    disp(subjInd)
    for sessInd = 1:2
        disp(subjInd)
        sess = ['sess' num2str(sessInd)];
        
%         x = main.dAll{subjInd}.(sess).data;
%         y = branch.dAll{subjInd}.(sess).data;
        
%         x = main.dAll{subjInd}.(sess).hr;
%         y = branch.dAll{subjInd}.(sess).hr;

        x = main.dAll{subjInd}.(sess).hr;
        y = branch.dAll{subjInd}.(sess).hr;

        
        ansTmp = ndims(x)==ndims(y)&&...
            all( size(x)==size(y) )&&...
            all( x(:)==y(:) );
        disp(ansTmp)
        ans = [ans ansTmp];
    end
end
all(ans)

