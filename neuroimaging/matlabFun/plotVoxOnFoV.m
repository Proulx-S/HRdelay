function plotVoxOnFoV(d,p,ind,flatFlag)
if ~exist('flatFlag','var')
    flatFlag = 0;
end
vec = mean(d.sin(ind,:),2);
[vec,~] = polarSpaceNormalization(vec,'cartRoi');
%         plotDensity(vec)
vec = abs(angle(vec));

cMap = redblue(256);
vecSpace = linspace(pi,0,256);
c = interp1(vecSpace,cMap,vec,'nearest');

figure('WindowStyle','docked');
if ~flatFlag
    ecc = d.voxProp.ecc(ind);
else
    ecc = d.voxProp.eccFlat(ind);
end
pol = d.voxProp.pol(ind)./180*pi;
hemiL = d.voxProp.hemifieldL(ind);
hemiR = d.voxProp.hemifieldR(ind);
% if flatFlag
%     polarscatter(pol(hemiL),log(ecc(hemiL)+1),eps,c(hemiL,:),'.','MarkerFaceColor','flat'); hold on
%     polarscatter(-pol(hemiR),log(ecc(hemiR)+1),eps,c(hemiR,:),'.');
% else
    polarscatter(pol(hemiL),ecc(hemiL),eps,c(hemiL,:),'.','MarkerFaceColor','flat'); hold on
    polarscatter(-pol(hemiR),ecc(hemiR),eps,c(hemiR,:),'.');
% end
ax = gca;
ax.ThetaZeroLocation = 'top';
% if flatFlag
%     ax.RTick = log((1:max(ecc))+1);
%     ax.RTickLabel = num2str((1:max(ecc))');
%     ax.RLim = [0 log(max(ecc)+1)];
% end

if ~isempty(p.featSel.fov.threshVal) && ~flatFlag
    th = repmat(linspace(0,2*pi,50),[2 1])';
%     if flatFlag
%         r = log(p.featSel.fov.threshVal+1);
%     else
        r = p.featSel.fov.threshVal;
%     end
    polarplot(th,r+zeros(size(th)),'k');
end


