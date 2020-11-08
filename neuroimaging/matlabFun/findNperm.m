function nPerm = findNperm(hitRate)

tmp = ~isnan(hitRate(:,any(~isnan(hitRate),1)));
if isempty(tmp)
    nPerm = 0;
else
    for i = 1:size(tmp,2);
        nPerm(i) = length(find(tmp(:,i)));
    end
end