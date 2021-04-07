function [featSel,f] = getFeatSel_areaAndFov(d,p)
warning('off','MATLAB:polyshape:repairedBySimplify')
visibleFlag = 0;

allFeatVal = cell(0);
allFeatP = cell(0);
% allFeatPrevInd = cell(0);
allFeatMethod = cell(0);
allFeatIndIn = cell(0);
condIndPairList = {[1 2 3]};

%% Precompute stats on response vector (random effect)
ind = true(size(d.sin,1),1);


%% Voxels representing stimulus fov
if p.featSel.fov.doIt
    curInfo1 = {'fov'};
    minContPercentArea = p.featSel.fov.empirical.minContPercentArea;
    thresh = p.featSel.(char(curInfo1));
    curInfo2 = {thresh.threshMethod};
    pVal = nan(size(ind));
    switch thresh.threshMethod
        case 'ecc'
            featVal = d.voxProp.ecc;
            curThresh = thresh.threshVal;
            
            curIndIn = curThresh(1)<featVal & featVal<curThresh(2);
        case 'empirical'
%             curFeatIndIn = allFeatIndIn(~cellfun('isempty',strfind(allFeatMethod,'vein:')) ...
%                 | ~cellfun('isempty',strfind(allFeatMethod,'act:')) ...
%                 | ~cellfun('isempty',strfind(allFeatMethod,'respVecSig:')));
%             prevInd = all(catcell(3,curFeatIndIn),3);
%             prevInd = prevInd(:,1);
            
            
            M = contourc(d.featSel.cont.X(1,:),d.featSel.cont.Y(:,1),double(d.featSel.cont.outXY),ones(2,1).*0.5);
            pgon = polyshape;
            while ~isempty(M)
                pgon = addboundary(pgon,M(1,2:1+M(2,1)),M(2,2:1+M(2,1)));
                M(:,1:1+M(2,1)) = [];
            end
            pgonOrig = pgon;
            
            pgon = regions(pgon);
            areas = area(pgon);
            pgon = pgon(areas/sum(areas) > minContPercentArea);
            
            X = d.featSel.cont.X;
            Y = d.featSel.cont.Y;
            U = d.featSel.cont.U;
            V = d.featSel.cont.V;
            tmpInd = false([size(X) length(pgon)]);
            tmpInd2 = false([size(U,1) length(pgon)]);
            for i = 1:length(pgon)
                Vertices = pgon(i).Vertices;
                tmpInd(:,:,i) = inpolygon(X,Y,Vertices(:,1),Vertices(:,2));
                tmpInd2(:,i) = inpolygon(U,V,Vertices(:,1),Vertices(:,2));
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
            imagesc(d.featSel.cont.X(1,:),d.featSel.cont.Y(:,1),~tmpInd); hold on
            set(gca,'YDir','normal'); colormap autumn
            plot(pgonOrig,'FaceColor','none')
            set(gca,'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1]);
            
            d.featSel.cont.outXY = ~tmpInd;
            prevInd = ind&tmpInd2;
            
            sm           = p.featSel.fov.empirical.auto(1).smList;
            mergeRadius  = p.featSel.fov.empirical.auto(1).mergeRadiusList;
            marginRadius = p.featSel.fov.empirical.auto(1).marginRadiusList;
            d = getDelayFovContour2(d,sm,prevInd);
            [~,f1,pgon] = processDelayFovContour2(d,p,prevInd,sm,mergeRadius,marginRadius,[],'do not add pgonRef',visibleFlag2);
            if p.featSel.fov.empirical.auto(2).smList~=sm
                sm           = p.featSel.fov.empirical.auto(2).smList;
                d = getDelayFovContour2(d,sm,prevInd);
            end
            mergeRadius  = p.featSel.fov.empirical.auto(2).mergeRadiusList;
            marginRadius = p.featSel.fov.empirical.auto(2).marginRadiusList;
            [curIndIn,f2,~] = processDelayFovContour2(d,p,prevInd,sm,mergeRadius,marginRadius,pgon,'add pgonRef',visibleFlag2);
            
            featVal = d.featSel.cont.vecUV;
        otherwise
            error('X')
    end
    f = [f0 f1 f2];
    curIndIn = repmat(curIndIn,[1 length(condIndPairList)]);
    
    allFeatVal(end+1) = {featVal};
    allFeatP(end+1) = {pVal};
    allFeatMethod(end+1) = {strjoin([curInfo1 curInfo2],': ')};
    allFeatIndIn(end+1) = {curIndIn};
    
%     ind = ind & curIndIn;
end


%% Multivariate (combined) feature selection
for i = 1:size(allFeatVal,2)
    allFeatVal{i} = permute(allFeatVal{i},[1 3 2]);
    allFeatP{i} = permute(allFeatP{i},[1 3 2]);
%     allFeatPrevInd{i} = permute(allFeatPrevInd{i},[1 3 2]);
    allFeatIndIn{i} = permute(allFeatIndIn{i},[1 3 2]);
end
sz = size(allFeatVal{1},[1 2 3]);
for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)>sz(3)
        sz(3) = size(allFeatVal{i},3);
    end
end
sz(2) = size(allFeatVal,2);

for i = 1:size(allFeatVal,2)
    if size(allFeatVal{i},3)==1
        allFeatVal{i} = repmat(allFeatVal{i},[1 1 sz(3)]);
    end
    if size(allFeatP{i},3)==1
        allFeatP{i} = repmat(allFeatP{i},[1 1 sz(3)]);
    end
%     if size(allFeatPrevInd{i},3)==1
%         allFeatPrevInd{i} = repmat(allFeatPrevInd{i},[1 1 sz(3)]);
%     end
    if size(allFeatIndIn{i},3)==1
        allFeatIndIn{i} = repmat(allFeatIndIn{i},[1 1 sz(3)]);
    end
end
featSel.featSeq.featVal = catcell(2,allFeatVal);
featSel.featSeq.featP = catcell(2,allFeatP);
% featSel.featSeq.featPrevInd = catcell(2,allFeatPrevInd);
featSel.featSeq.featQtile = nan(size(featSel.featSeq.featVal));
featSel.featSeq.featIndIn = catcell(2,allFeatIndIn);
featSel.featSeq.featSelList = allFeatMethod;
featSel.featSeq.condPairList = permute(condIndPairList,[1 3 2]);
featSel.featSeq.info = 'vox X featSel X condPair';


switch p.featSel.global.method
    case 'all'
        indIn = all(featSel.featSeq.featIndIn,2);
        % Compute quantile
        for featInd = 1:sz(2)
            for condIndPairInd = 1:sz(3)
                x = featSel.featSeq.featVal(:,featInd,condIndPairInd);
                prevInd = true(size(x));
                [fx,x2] = ecdf(x(prevInd));
                x2 = x2(2:end); fx = fx(2:end);
                [~,b] = ismember(x,x2);
                featSel.featSeq.featQtile(b~=0,featInd,condIndPairInd) = fx(b(b~=0));
            end
        end
    otherwise
        error('code that')
        
        allFeatVal = catcell(2,allFeatVal);
        x = zscore(-log(allFeatVal(:,1)));
        y = zscore(log(allFeatVal(:,2)));
        z = zscore(log(allFeatVal(:,3)));
        
        
        [fx,x2] = ecdf(x);
        x2 = x2(2:end); fx = fx(2:end);
        [~,b] = ismember(x,x2);
        fx = fx(b); clear x2
        
        [fy,y2] = ecdf(y);
        y2 = y2(2:end); fy = fy(2:end);
        [~,b] = ismember(y,y2);
        fy = fy(b); clear y2
        
        [fz,z2] = ecdf(z);
        z2 = z2(2:end); fz = fz(2:end);
        [~,b] = ismember(z,z2);
        fz = fz(b); clear z2
        
        %     fyz = fy.*fz;
        %     fyz = fyz - min(fyz);
        %     fyz = fyz./max(fyz);
        fxyz = fx.*fy.*fz;
        fxyz = fxyz - min(fxyz);
        fxyz = fxyz./max(fxyz);
        
        perc = p.featSel.global.percentile;
        %     ind2 = fyz>prctile(fyz,perc) & fx>prctile(fx,perc);
        ind2 = fxyz>prctile(fxyz,perc) & fx>prctile(fx,perc);
        
        %     if p.figOption.verbose>=1 && subjInd==p.figOption.subjInd && sessInd==p.figOption.sessInd
        %         figure('WindowStyle','docked');
        %         scatter3(x(ind2),y(ind2),z(ind2),'k.'); hold on
        %         scatter3(x(~ind2),y(~ind2),z(~ind2),'r.');
        %         % ax = gca;
        %         % ax.CameraPosition = camPos;
        %         xlabel(['zscore(-log(' info1{1+1} '))'])
        %         ylabel(['zscore(log(' info1{1+2} '))'])
        %         zlabel(['zscore(log(' info1{1+3} '))'])
        %         legend({'include' 'exclude'})
        %         % ax = gca; camPos = ax.CameraPosition;
        %         % ax = gca; ax.CameraPosition = camPos;
        %     end
end
%% Output
featSel.indIn = indIn;
featSel.info = 'vox X featSel X condPair';

