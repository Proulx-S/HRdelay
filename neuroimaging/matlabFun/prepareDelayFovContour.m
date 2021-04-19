function cont = prepareDelayFovContour(d,voxProp,p)
padFac = p.featSel.fov.empirical.padFac;
cont = cell(size(d));
for subjInd = 1:size(d,1)
    curVoxProp = voxProp{subjInd};
    % Apply precomputed flattening
    if isfield(curVoxProp,'eccTrans')
        curVoxProp.ecc = curVoxProp.eccTrans(curVoxProp.ecc);
    end
    for sessInd = 1:size(d,2)
        vecUV = mean(d{subjInd,sessInd}.sin(curVoxProp.hemifield,:),2);
        
        % Prepare surface
        [~,U,V,densityXY,X,Y] = pol2surf(curVoxProp,padFac);
        outXY = isnan(densityXY); clear densityXY
        [vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
        vecUV = abs(angle(vecUV));
        
        % Output
        cont{subjInd,sessInd}.U = U;
        cont{subjInd,sessInd}.V = V;
        cont{subjInd,sessInd}.vecUV = vecUV;
        cont{subjInd,sessInd}.X = X;
        cont{subjInd,sessInd}.Y = Y;
        cont{subjInd,sessInd}.outXY = outXY;
    end
end
