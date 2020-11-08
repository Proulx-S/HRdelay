function [circT,alphaDelta] = getCircT(data)



alpha = sec2rad(data,12);
alpha = [alpha(:,:,1) alpha(:,:,2)];
% alpha = [alpha(:,randperm(size(alpha,2)),1) alpha(:,randperm(size(alpha,2)),2)];
idx = [ones(size(data,2),1)*1; ones(size(data,2),1)*2];

for i = 1:size(alpha,1)
    alpha1 = alpha(i,idx==1);
    alpha2 = alpha(i,idx==2);
    alphaDelta(i) = circ_mean(alpha1')-circ_mean(alpha2');
    [~, table{i}] = circ_wwtest(alpha1,alpha2);
    F(i) = table{i}{2,5};
    MSerror(i) = table{i}{3,4};
    circT(i) = alphaDelta(i)/sqrt(2*MSerror(i)/length(alpha1));
end
