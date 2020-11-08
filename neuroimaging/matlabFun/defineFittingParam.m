function p = defineFittingParam(ts,ind,data,fixOnly)
if ~exist('fixOnly','var')
    fixOnly = 0;
end

p.L = size(ts,2);
p.deltaT = 1;
p.period = 12;
p.crossCorrGrain = 0.05;
p.dataName = data;
p.ind = ind;
p.fixDelay = 1;

%Define time
p.t = 0:p.deltaT:(p.L-1)*p.deltaT;
p.tShifts = 0:p.crossCorrGrain:p.period/2-p.crossCorrGrain;
for i = 1:length(p.tShifts)
    p.model(i,:) = -cos(2*pi*1/p.period*(p.t-p.tShifts(i)));
end

if fixOnly
    [~,tmpInd] = min(abs(p.tShifts-p.fixDelay));
    p.model = p.model(tmpInd,:);
    p.tShifts = p.tShifts(tmpInd);
end