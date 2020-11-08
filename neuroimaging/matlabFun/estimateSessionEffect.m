function [betaHat,X] = estimateSessionEffect(data,sessionLabel)
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
betaHat = nan(size(data,2),size(X,2));
for vox = 1:size(data,2)
    A = data(:,vox,1);
    B = data(:,vox,2);
    C = data(:,vox,3);
    %contruct data matrix y
    y = cat(1,A,B,C);
    %solve with OLS
    betaHat(vox,:) = (X' * X) \ X' * y; % computing projection of matrix X on y, giving beta
%     betaHat = inv(X' * X) * X' * y; % computing projection of matrix X on
%     y, giving beta (should be the same as above)
%     betaHat = pinv(X' * X) * X' * y; % use pinv because singular to machine precision
end



