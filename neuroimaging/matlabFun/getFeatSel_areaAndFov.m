function [featVal,featMethod,featIndIn,featInfo,f] = getFeatSel_areaAndFov(cont,voxProp,p)

warning('off','MATLAB:polyshape:repairedBySimplify')
visibleFlag = 0;

%% Precompute stats on response vector (random effect)
ind = true(size(cont.U,1),1);

%% Voxels representing stimulus fov
if p.featSel.fov.doIt
    curInfo1 = {'fov'};
    minContPercentArea = p.featSel.fov.empirical.minContPercentArea;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
%     pVal = nan(size(ind));
    switch thresh.threshMethod
        case 'ecc'
            error('double-check that')
            featVal = d.voxProp.ecc;
            curThresh = thresh.threshVal;
            
            curIndIn = curThresh(1)<featVal & featVal<curThresh(2);
        case 'empirical'
            M = contourc(cont.X(1,:),cont.Y(:,1),double(cont.outXY),ones(2,1).*0.5);
            pgon = polyshape;
            while ~isempty(M)
                pgon = addboundary(pgon,M(1,2:1+M(2,1)),M(2,2:1+M(2,1)));
                M(:,1:1+M(2,1)) = [];
            end
            pgonOrig = pgon;
            
            pgon = regions(pgon);
            areas = area(pgon);
            pgon = pgon(areas/sum(areas) > minContPercentArea);
            
            tmpInd = false([size(cont.X) length(pgon)]);
            tmpInd2 = false([size(cont.U,1) length(pgon)]);
            for i = 1:length(pgon)
                Vertices = pgon(i).Vertices;
                tmpInd(:,:,i) = inpolygon(cont.X,cont.Y,Vertices(:,1),Vertices(:,2));
                tmpInd2(:,i) = inpolygon(cont.U,cont.V,Vertices(:,1),Vertices(:,2));
            end
            tmpInd = any(tmpInd,3);
            tmpInd2 = any(tmpInd2,2);
            if visibleFlag
                visibleFlag2 = 2;
            else
                visibleFlag2 = 1;
            end
            if visibleFlag
                f0 = figure('WindowStyle','docked','visible','on');
            else
                f0 = figure('WindowStyle','docked','visible','off');
            end
            imagesc(cont.X(1,:),cont.Y(:,1),~tmpInd); hold on
            set(gca,'YDir','normal'); colormap autumn
            plot(pgonOrig,'FaceColor','none')
            set(gca,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1]);
            
            cont.outXY = ~tmpInd;
            prevInd = ind&tmpInd2;
            
            sm           = p.featSel.fov.empirical.auto(1).smList;
            mergeRadius  = p.featSel.fov.empirical.auto(1).mergeRadiusList;
            marginRadius = p.featSel.fov.empirical.auto(1).marginRadiusList;
            cont = getDelayFovContour3(cont,sm,prevInd);
            featVal = cont.vecUV;
            [~,f1,pgon] = processDelayFovContour3(cont,voxProp,p,prevInd,sm,mergeRadius,marginRadius,[],'do not add pgonRef',visibleFlag2);
            
            if p.featSel.fov.empirical.auto(2).smList~=sm
                sm = p.featSel.fov.empirical.auto(2).smList;
                cont = getDelayFovContour3(cont,sm,prevInd);
            end
            mergeRadius  = p.featSel.fov.empirical.auto(2).mergeRadiusList;
            marginRadius = p.featSel.fov.empirical.auto(2).marginRadiusList;
            if ~isempty(pgon)
                [curIndIn,f2,~] = processDelayFovContour3(cont,voxProp,p,prevInd,sm,mergeRadius,marginRadius,pgon,'add pgonRef',visibleFlag2);
            else
                curIndIn = true(size(prevInd));
                f2 = gobjects([1 5]);
            end
        otherwise
            error('X')
    end
    f = [f0 f1 f2];
end
featVal;
featMethod = strjoin([curInfo1 curInfo2],': ');
featIndIn = curIndIn;
featInfo = 'vox X featSel X condPair';
