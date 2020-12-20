clear all
close all

main = load('C:\Users\sebas\Desktop\main.mat');
branch = load('C:\Users\sebas\Desktop\branch.mat');

ans = [];
for subjInd = 1:length(main.d)
    for sessInd = 1:2
        disp(['subj' num2str(subjInd) ' sess' num2str(sessInd)])
        sess = ['sess' num2str(sessInd)];
        
%         x = main.d{subjInd}.(sess).data;
%         y = branch.dC{subjInd}.(sess).data;
        
%         x = main.d{subjInd}.(sess).hr;
%         y = branch.dC{subjInd}.(sess).hr;
        
%         x = main.d{subjInd}.(sess).runLabel;
%         y = branch.dC{subjInd}.(sess).runLabel;
        
%         x = main.featSelStats{subjInd}.anyCondActivation.(sess).F;
%         y = branch.dC{subjInd}.(sess).anyCondActivation_F;
        
%         x = main.featSelStats{subjInd}.anyCondActivation.(sess).FDR;
%         y = branch.dC{subjInd}.(sess).anyCondActivation_FDR;
        
%         x = main.featSelStats{subjInd}.anyCondActivation.(sess).P;
%         y = branch.dC{subjInd}.(sess).anyCondActivation_P;
        
        x = main.vein{subjInd}.(sess).noiseOverMean;
        y = branch.dC{subjInd}.(sess).vein_score;
        
        x = main.vein{subjInd}.sess1.noiseOverMean;
        y = branch.dC{subjInd}.sess2.vein_score;
        
        
%         x = main.vein{subjInd}.(sess).thresh;
%         y = branch.dC{subjInd}.(sess).vein_thresh;

        scatter(x,y)
        xlabel('main branch vein score')
        ylabel('move-all... branch vein score')
        
        
%         if subjInd==2 && sessInd==1
%             x(4,:,:,:) = [];
%         end
        ansTmp = ndims(x)==ndims(y)&&...
            all( size(x)==size(y) )&&...
            all( x(:)==y(:) );
        disp(ansTmp)
        ans = [ans ansTmp];
    end
end
all(ans)



clear all
close all

load('C:\Users\sebas\Desktop\branchVein.mat');
sum(runInd & sessLabel==sessInd)
featSel.(sess).vein_score(maskFitArea) = mean(results.OLS.mixed.veinFull(:,:,:,runInd & sessLabel==sessInd),4);
imagesc(featSel.(sess).vein_score(:,:,10))

load('C:\Users\sebas\Desktop\mainVein.mat');
vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFitArea) = mean(results.OLS.mixed.veinFull(:,:,:,ind),4);
