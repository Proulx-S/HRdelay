function [x,y,k] = getXYK(d,p)
if ~isfield(p,'condPair')
    p.condPair = 'all';
end
switch p.condPair
    case 'grat1VSgrat2'
        condInd = [1 2];
    case 'grat1VSplaid'
        condInd = [1 3];
    case 'grat2VSplaid'
        condInd = [2 3];
    case 'all'
        condInd = [1 2 3];
    otherwise
        error('code that')
end
% Define x(data), y(label) and k(xValFolds)
x = cell(size(condInd));
y = cell(size(condInd));
k = cell(size(condInd));
nRep = size(d.sin(:,:,:,:,condInd,:),4);
for i = 1:length(condInd)
    x{i} = permute(d.sin(:,:,:,:,condInd(i),:),[4 1 2 3]);
end
for i = 1:length(condInd)
    y{i} = i.*ones(nRep,1);
    k{i} = (1:nRep)';
end
x = catcell(1,x);
y = catcell(1,y);
k = catcell(1,k);