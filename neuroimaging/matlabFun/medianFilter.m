function [xbin,ymed] = medianFilter(x,y,n)

[~,edges,bin] = histcounts(x,n);
xbin = edges(1:end-1)+diff(edges([1 end]))/n/2;
ymed = nan(size(xbin));
for binInd = 1:n
    ymed(binInd) = median(y(bin==binInd));
end
