function [densityUV,U,V,densityXY,X,Y] = pol2surf(voxProp,padFac,n)
if ~exist('padFac','var') || isempty(padFac)
    padFac = 1.2;
end
if ~exist('n','var') || isempty(n)
    n = 2^8;
end
ecc = voxProp.ecc;
pol = voxProp.pol;
pol(voxProp.hemifieldL) = voxProp.pol(voxProp.hemifieldL)./180*pi;
pol(voxProp.hemifieldR) = -voxProp.pol(voxProp.hemifieldR)./180*pi;
pol = wrapToPi(pol+pi/2);
[U,V] = pol2cart(pol,ecc);

bwFac = 0.009;
delta = range([U; V])*padFac/n;
[X,Y] = meshgrid(linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),n);
[densityXY,~] = ksdensity([U V],[X(:) Y(:)],'Bandwidth',delta*n*bwFac);
densityXY = reshape(densityXY,[n n]);
densityUV = interp2(X,Y,densityXY,U,V);
densityXY(densityXY<(min(densityUV)/10)) = nan;


