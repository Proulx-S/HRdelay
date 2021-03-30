function processFeatSel(p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            outDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,2)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
%     for sessInd = 1:2
%         sess = ['sess' num2str(sessInd)];
%         if isfield(res.(sess),'featSel')
%             res.(sess) = rmfield(res.(sess),'featSel');
%         end
%     end
    dAll{subjInd} = res;
end
d = dAll; clear dAll
sessList = fields(d{1});


%% Reorganize
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        sess = ['sess' num2str(sessInd)];
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
    end
end
d = dP; clear dP

%% Empirical FOV
% p.figOption.verbose
d = flattenEccDist(d,p,1);

for subjInd = 2:size(d,1)
    sessInd = 1;
    
    getDelayFovMap(d{subjInd,sessInd},1);
    
    sessInd = 1;
    voxProp = d{subjInd,sessInd}.voxProp;
    voxProp.ecc = voxProp.eccTrans(voxProp.ecc);
    [~,U,V,densityXY,X,Y] = pol2surf(voxProp);
    
    ind = true(size(d{subjInd,sessInd}.sin,1),1);
    vecUV = mean(d{subjInd,sessInd}.sin(ind,:),2);
    [vecUV,~] = polarSpaceNormalization(vecUV,'cartRoi');
    vecUV = abs(angle(vecUV));
    F = scatteredInterpolant(U,V,vecUV,'natural','nearest');
    vecXY = F(X,Y);
%     vecXY(isnan(densityXY)) = nan;
    cMap = redblue(256);
    vecSpace = linspace(pi,0,256);
    C = interp1(vecSpace,cMap,vecXY,'nearest');
    
    
    
    
    
    figure('WindowStyle','docked');
%     imagesc(X(1,:),Y(:,1),imgaussfilt(vecXY),[0 pi]); hold on
    imagesc(X(1,:),Y(:,1),vecXY,[0 pi]); hold on
    colormap(flip(cMap,1));
    set(gca,'YDir','normal')
    hScat = scatter(U,V);
    hScat.MarkerEdgeColor = 'none';
    hScat.MarkerFaceColor = 'k';
    hScat.SizeData = 4^2;
    hScat.MarkerFaceAlpha = 0.1;
    ax = gca;
    ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
    ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
    xLim = xlim; yLim = ylim;
    
    [M,c] = contour(X,Y,imgaussfilt(vecXY,1),[1 1].*pi/2); hold on
    c.LineColor = 'k';
%     c.ZData(:) = max(vecXY(:))+1;
%     scatter3(M(1,2:end),M(2,2:end),ones(1,size(M(:,2:end),2)).*max(vecXY(:))+1)
%     plot3(M(1,2:end),M(2,2:end),ones(1,size(M(:,2:end),2)).*max(vecXY(:))+1)
    
    
    
    
    figure('WindowStyle','docked');
    surf(X,Y,vecXY,C,'LineStyle','none'); hold on
    view([0,90]);
    scatter3(U,V,vecUV,eps,'k.')
    ylim(yLim); xlim(xLim);
    ax = gca;
    ax.DataAspectRatio = ax.DataAspectRatio([1 1 3]);
    ax.PlotBoxAspectRatio = ax.PlotBoxAspectRatio([1 1 3]);
    
    [M,c] = contour(X,Y,vecXY,[1 1].*pi/2); hold on
    c.ZData(:) = max(vecXY(:))+1;
    scatter3(M(1,2:end),M(2,2:end),ones(1,size(M(:,2:end),2)).*max(vecXY(:))+1)
    plot3(M(1,2:end),M(2,2:end),ones(1,size(M(:,2:end),2)).*max(vecXY(:))+1)
    
    
    
    for i = 1:2
        ecc = voxProp.eccTrans(p.featSel.fov.threshVal(i));
        pol = linspace(-pi,pi,100);
        [u,v] = pol2cart(pol,ecc);
        plot3(u,v,max(vecUV).*ones(size(v)),'k')
    end
    eccList = 0:6;
    for eccInd = 1:length(eccList)
        ecc = voxProp.eccTrans(eccList(eccInd));
        pol = linspace(-pi,pi,100);
        [u,v] = pol2cart(pol,ecc);
        plot3(u,v,max(vecUV).*ones(size(v)),'color',[1 1 1].*0.8)
    end
    
    
end


%% Feature selection
featSel = cell(size(d));
disp('computing feature selection stats')
for sessInd = 1:numel(d)
    disp(['sess' num2str(sessInd) '/' num2str(numel(d))])
    [subjInd,sessInd] = ind2sub(size(d),sessInd);
    featSel{subjInd,sessInd} = getFeatSel(d{subjInd,sessInd},p);
end
disp('done')


indInX = cell(1,size(d,2));
dX = cell(1,size(d,2));
fieldList = fields(d{subjInd,sessInd});
for sessInd = 1:size(d,2)
    for subjInd = 1:size(d,1)
        if subjInd==1
            indInX{sessInd} = featSel{subjInd,sessInd}.featSeq.featIndIn;
            dX{sessInd} = d{subjInd,sessInd};
        else
            indInX{sessInd} = cat(1,indInX{sessInd},featSel{subjInd,sessInd}.featSeq.featIndIn);
            for fieldInd = 1:length(fieldList)
                if isnumeric(d{subjInd,sessInd}.(fieldList{fieldInd}))...
                        && ~strcmp(fieldList{fieldInd},'sinDesign')...
                        && ~strcmp(fieldList{fieldInd},'hrDesign')
                    dX{sessInd}.(fieldList{fieldInd}) = cat(1,mean(dX{sessInd}.(fieldList{fieldInd}),4),mean(d{subjInd,sessInd}.(fieldList{fieldInd}),4));
                elseif isstruct(d{subjInd,sessInd}.(fieldList{fieldInd}))...
                        && ~strcmp(fieldList{fieldInd},'featSel')
                    fieldList2 = fields(d{subjInd,sessInd}.(fieldList{fieldInd}));
                    for fieldInd2 = 1:length(fieldList2)
                        if isnumeric(d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}))...
                                || islogical(d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}))
                            dX{sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}) = ...
                                cat(1,dX{sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}),d{subjInd,sessInd}.(fieldList{fieldInd}).(fieldList2{fieldInd2}));    
                        end
                    end
                end
            end
        end
    end
end

sessInd = 1;
plotVoxOnFoV(dX{sessInd},p,true(size(dX{sessInd}.sin,1),1))
ax = gca;
RLim = ax.RLim;
title('allVox')
for featInd = 1:length(featSel{subjInd,sessInd}.featSeq.featSelList)
    plotVoxOnFoV(dX{sessInd},p,indInX{sessInd}(:,featInd,1))
    title(featSel{subjInd,sessInd}.featSeq.featSelList{featInd})
    ax = gca;
    ax.RLim = RLim;
end

featSel{subjInd,sessInd}.featSeq.featSelList'
featInd = [3 4];
plotVoxOnFoV(dX{sessInd},p,all(indInX{sessInd}(:,featInd,1),2))
ax = gca;
R = [ax.Children(3).RData ax.Children(4).RData];
figure('WindowStyle','docked');
R = exp(R)-1;
hist(R,100)




featInd = [0];
condPairInd = 1;
if featInd==0
    ind = true(size(featSel{subjInd,sessInd}.featSeq.featIndIn,1),1);
else
    ind = all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,featInd,condPairInd),2);
end
p.featSel.fov.threshVal = [];
plotVoxOnFoV(dX{1},p,ind)

% fac = 1;
% fac = (2*pi);
% fac = 1/(2*pi);
% fac = (2*2*pi);
fac = 1/(2*2*pi);

R = dX{sessInd}.voxProp.ecc;
Rp = cdf(nonparamDistFit(R),R);
v = linspace(min(R),max(R)*fac,length(unique(R)));%R
x = linspace(1/length(unique(R)),1,length(unique(R)));%Rp
xq = Rp;%Rp
vq = interp1(x,v,xq);%R
R2 = vq;
Rp2 = xq;
figure('WindowStyle','docked');
[~,uInd,~] = unique(R);
plot(R(uInd),Rp(uInd),'.'); hold on
[~,uInd,~] = unique(R2);
plot(R2(uInd),Rp2(uInd),'.'); hold on
nonparamDistFit(R2,0)

d2 = dX{sessInd};
d2.voxProp.ecc = R2;
plotVoxOnFoV(d2,p,ind)





%% Save
disp('saving feature selection')
if ~exist(fullfile(funPath,outDir),'dir')
    mkdir(fullfile(funPath,outDir))
end
fullfilename = fullfile(funPath,outDir,'featSel.mat');
save(fullfilename,'featSel')
disp(['saved to: ' fullfilename])
