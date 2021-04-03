function d = prepareDelayFovContour(d,p)
padFac = p.featSel.fov.empirical.padFac;
for sessInd = 1:numel(d)
    voxProp = d{sessInd}.voxProp;
    vecUV = mean(d{sessInd}.sin(:,:),2);
    
    % Apply precomputed flattening
    if isfield(voxProp,'eccTrans')
        voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
    end
    
    % Prepare surface
    [~,U,V,densityXY,X,Y] = pol2surf(voxProp,padFac);
    outXY = isnan(densityXY); clear densityXY
    [vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
    vecUV = abs(angle(vecUV));
    
    % Output
    d{sessInd}.featSel.cont.U = U;
    d{sessInd}.featSel.cont.V = V;
    d{sessInd}.featSel.cont.vecUV = vecUV;
    d{sessInd}.featSel.cont.X = X;
    d{sessInd}.featSel.cont.Y = Y;
    d{sessInd}.featSel.cont.outXY = outXY;
end
