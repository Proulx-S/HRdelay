clear all
close all

main = load('C:\Users\sebas\Desktop\main.mat');
branch = load('C:\Users\sebas\Desktop\branch.mat');


for subjInd = 1:length(main)
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        main.dC{subjInd}.(sess).data==branch.dC{subjInd}.(sess).data
        
        
    end
end

