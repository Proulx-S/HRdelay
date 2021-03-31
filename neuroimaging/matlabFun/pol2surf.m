function [densityUV,U,V,densityXY,X,Y] = pol2surf(voxProp)
ecc = voxProp.ecc;
pol = voxProp.pol;
pol(voxProp.hemifieldL) = voxProp.pol(voxProp.hemifieldL)./180*pi;
pol(voxProp.hemifieldR) = -voxProp.pol(voxProp.hemifieldR)./180*pi;
pol = wrapToPi(pol+pi/2);
[U,V] = pol2cart(pol,ecc);

n = 2^7;
bwFac = 0.009;
delta = range([U; V])*1.2/n;
[X,Y] = meshgrid(linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),linspace(min([U; V])-0.1*delta*n,max([U; V])+0.1*delta*n,n),n);
[densityXY,~] = ksdensity([U V],[X(:) Y(:)],'Bandwidth',delta*n*bwFac);
densityXY = reshape(densityXY,[n n]);
densityUV = interp2(X,Y,densityXY,U,V);
densityXY(densityXY<(min(densityUV)/10)) = nan;


