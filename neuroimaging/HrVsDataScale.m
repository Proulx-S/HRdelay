clear all
close all
load tmp

hr = squeeze(mean(tmpHr,3))'; clear tmpHr
data = squeeze(mean(tmpData,3)); clear tmpData


t = 1:12;
for voxInd = 1:size(data,2)
    [fitresult, ~] = sinFit(t, hr(:,voxInd)');
    a(voxInd) = fitresult.a;
    c(voxInd) = fitresult.c;
end

figure('WindowStyle','docked');
scatter(a(a>0),abs(data(a>0)))
lsline
polyfit(a(a>0),abs(data(a>0)),1)
xlabel('hr')
ylabel('data')



