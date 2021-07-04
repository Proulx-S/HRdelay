function [eccL,eccR,densityL,densityR] = compareEcc(d,eccTransL,eccTransR)
dOrig = d;
sm = 0.99; plotFlag = 0;

%% Left
d = dOrig;
d.voxProp.ecc = d.voxProp.ecc(d.voxProp.hemifieldL);
d.voxProp.pol = d.voxProp.pol(d.voxProp.hemifieldL);
d.voxProp.hemifieldR = d.voxProp.hemifieldR(d.voxProp.hemifieldL);
d.voxProp.hemifieldL = d.voxProp.hemifieldL(d.voxProp.hemifieldL);

ecc = d.voxProp.ecc;
if exist('eccTransL','var') && ~isempty(eccTransL)
    ecc = eccTransL.toFlat{end}(ecc);
end
[~,~,~,~,~,~,~,eccXYspline,densityXYspline,~,~] = getEccPd(d,ecc,sm,plotFlag);
if plotFlag
    figure('WindowStyle','docked');
    plot(eccXYspline,densityXYspline)
end
eccL = eccXYspline;
densityL = densityXYspline;




%% Left
d = dOrig;
d.voxProp.ecc = d.voxProp.ecc(d.voxProp.hemifieldR);
d.voxProp.pol = d.voxProp.pol(d.voxProp.hemifieldR);
d.voxProp.hemifieldL = d.voxProp.hemifieldL(d.voxProp.hemifieldR);
d.voxProp.hemifieldR = d.voxProp.hemifieldR(d.voxProp.hemifieldR);

ecc = d.voxProp.ecc;
if exist('eccTransR','var') && ~isempty(eccTransR)
    ecc = eccTransR.toFlat{end}(ecc);
end
[~,~,~,~,~,~,~,eccXYspline,densityXYspline,~,~] = getEccPd(d,ecc,sm,plotFlag);
if plotFlag
    figure('WindowStyle','docked');
    plot(eccXYspline,densityXYspline)
end
eccR = eccXYspline;
densityR = densityXYspline;

