function cont = prepareDelayFovContour(d,voxProp,p)
padFac = p.featSel.fov.empirical.padFac;
cont = cell(size(d));
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        vecUV = mean(d{subjInd,sessInd}.sin(voxProp{subjInd}.hemifield,:),2); d{subjInd,sessInd} = {};
        
        % Apply precomputed flattening
        if isfield(voxProp{subjInd},'eccTrans')
            voxProp{subjInd}.ecc = voxProp{subjInd}.eccTrans(voxProp{subjInd}.ecc);
        end
        
        % Prepare surface
        [~,U,V,densityXY,X,Y] = pol2surf(voxProp{subjInd},padFac);
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
