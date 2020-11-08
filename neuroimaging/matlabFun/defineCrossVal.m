
function crossVal = defineCrossVal(d,p)

nSplit = str2double(p.split(1));
tmp1 = [];
tmp2 = [];
for i = 1:nSplit
    tmp1 = [tmp1 randperm(size(d.label,1)/2/nSplit)];
    tmp2 = [tmp2 randperm(size(d.label,1)/2/nSplit)];
end
tmp1 = tmp1(randperm(length(tmp1)));
tmp2 = tmp2(randperm(length(tmp2)));
crossVal = cat(1,tmp1',tmp2');


