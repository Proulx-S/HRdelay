function contData = processDelayFovContour(d,p,contData,contInd,subjInd,interp,extrap,mergeRadius,marginRadius)
warning('off','MATLAB:polyshape:repairedBySimplify')

F = scatteredInterpolant(contData.U,contData.V,contData.vecUV,interp,extrap);
vecXY = F(contData.X,contData.Y);
vecXY(contData.outXY|contData.outUV) = pi/2;
showDelayFovContour(contData.X,contData.Y,vecXY,contData.cMap,contData.U,contData.V,[],[],'on')
addEccRef(d,p)
title(['subj' num2str(subjInd)])
drawnow

% mergeRadius = 1;
% marginRadius = 0.5;
%% Process contours
if contInd==inf
    contInd = 1:length(contData.cont);
end
% Compile into polygons
pgon = polyshape;
for i = 1:length(contInd)
    x = contData.cont{contInd(i)}(:,1);
    y = contData.cont{contInd(i)}(:,2);
    pgon = addboundary(pgon,x,y);
end
% Remove holes (contour inside another contour)
pgon = rmholes(pgon);
plot(pgon,'FaceAlpha',0);
% Connect nearby contours by expanding then shrinking them
pgon = polybuffer(pgon,mergeRadius);
pgon = polybuffer(pgon,-mergeRadius);
pgon = simplify(pgon);
plot(pgon,'FaceAlpha',0);
% Expand for conservative boundaries
pgon = polybuffer(pgon,marginRadius);
plot(pgon,'FaceAlpha',0);

pgon = regions(pgon);
[~,b] = sort(area(pgon),'descend');
pgon = pgon(b);

contData.contPgon = polyshape;
for i = 1:length(pgon)
    contData.contPgon = addboundary(contData.contPgon,pgon(i).Vertices(:,1),pgon(i).Vertices(:,2));
end
    