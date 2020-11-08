function [dataClean,X,curX,allBeta] = regressOutSessionVoxWise2(data,sessionLabel,plotIt)
% grat1 = grat1' + X;
% grat2 = grat2' + X;
% plaid = plaid' + X;
% is an overdetermined system, but approximation can be found using
% ordinary least squares on the 11 measures, performed independently at
% each voxel
%https://en.wikipedia.org/wiki/Overdetermined_system#Approximate_solutions
%https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)#Derivation_of_the_normal_equations

% grat = abs(d.xData)';
% plaid = abs(d.normData)';
% grat_xcorr_before = corr(grat);
% grat_xcorr_before(logical(diag(ones(size(grat_xcorr_before,1),1)))) = nan;

% data = cat(3,grat(:,1:end/2),grat(:,end/2+1:end),plaid);
dataClean = nan(size(data));
%define design matrix X
XA = [1 0 1];
XB = [0 1 1];
XC = [0 0 1];
Xs = sessionLabel;
tmp = [];
for i = 1:length(unique(Xs))
    tmp = cat(2,tmp,cat(1,Xs(:,:,1),Xs(:,:,2),Xs(:,:,3))==i);
end
Xs = tmp(:,2:end); clear tmp

X1 = cat(1,repmat(XA,[size(data,1) 1]),repmat(XB,[size(data,1) 1]),repmat(XC,[size(data,1) 1]));
X = cat(2,X1,Xs); %such that columns represent the following regressors [grat1 grat2 plaid common session1 session2 session3]
allBeta = nan(size(data,2),size(X,2));
for vox = 1:size(data,2)
    A = data(:,vox,1);
    B = data(:,vox,2);
    C = data(:,vox,3);
    %contruct data matrix y
    y = cat(1,A,B,C);
    %solve with OLS
    betaHat = (X' * X) \ X' * y; % computing projection of matrix X on y, giving beta
%     betaHat = inv(X' * X) * X' * y; % computing projection of matrix X on
%     y, giving beta (should be the same as above)
%     betaHat = pinv(X' * X) * X' * y; % use pinv because singular to machine precision
    %recontruct the nuisance response
    curX = X; curX(:,1:3) = 0;
    nuisance = curX*betaHat;
    %subtract nuisance from data
    y_clean = y-nuisance;
    %recon data back to original structure
    dataClean(:,vox,:) = reshape(y_clean,[size(dataClean,1) 1 size(dataClean,3)]);
    allBeta(vox,:) = betaHat;
end



% 
% 
% grat = dataClean(:,:,[1 2]);
% grat = cat(2,grat(:,:,1),grat(:,:,2));
% grat_xcorr_after = corr(grat);
% grat_xcorr_after(logical(diag(ones(size(grat_xcorr_after,1),1)))) = nan;
% 
% 
% if plotIt
%     
%     figure('windowStyle','docked');
%     imagesc(grat_xcorr_before,[-1 1]); colorbar
%     saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '_1bf.jpg']))
%     figure('windowStyle','docked');
%     imagesc(grat_xcorr_after,[-1 1]); colorbar
%     saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '_1af.jpg']))
%     
%     lims = [min([grat_xcorr_before(:); grat_xcorr_after(:)]) max([grat_xcorr_before(:); grat_xcorr_after(:)])];
%     figure('windowStyle','docked');
%     imagesc(grat_xcorr_before,lims); colorbar
%     saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '_2bf.jpg']))
%     figure('windowStyle','docked');
%     imagesc(grat_xcorr_after,lims); colorbar
%     saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '_2af.jpg']))
% end
% 
% 
% %% Put back in appropriate output format
% dataAmp = cat(1,cat(2,dataClean(:,:,1),dataClean(:,:,2))',plaid');
% dataPhase = angle(cat(1,d.xData,d.normData));
% [X,Y] = pol2cart(dataPhase,dataAmp);
% dataClean = complex(X,Y);
% d.xData = dataClean(1:end*2/3,:);
% d.normData = dataClean(end*2/3+1:end,:);

