function [fFOV,fDensity] = replotOneFOV(p,funPath,inDir,subjInd,sessInd)


fFOV = openfig(fullfile(p.figOption.finalDir,'processFOV.fig'));
ax = findobj(fFOV.Children,'Type','axes');
tmp = cell(size(ax));
for i = 1:length(ax)
    tmp(i) = ax(i).Title.String(1);
end
fFOV2 = figure('WindowStyle','docked');
ax = copyobj(ax(ismember(tmp,['subj' num2str(subjInd) '; sess' num2str(sessInd)])),fFOV2);
drawnow
close(fFOV)
fFOV = fFOV2; clear fFOV2
fFOV.Units = 'inches';
set(ax,'Units','centimeter')
drawnow
ax(1).Position = [7    1.6845   10.1600   11.6152];
ax(2).Position = [0    1.6845   10.1600   11.6152];

delete(ax(1).Title)
delete(ax(2).Title)
ax(1).YAxis.Visible = 'on';
ax(2).YAxis.Visible = 'on';
ax(2).YAxisLocation = 'right';

% hScat = findobj(ax(1).Children,'Type','Scatter');
% hScat.MarkerEdgeColor
% hScat.MarkerFaceAlpha

fData = load([p.featSel.fov.resFile,'.mat']);

eccTrans = fData.voxProp{subjInd,sessInd}.L.eccTrans;
ecc = fData.voxProp{subjInd,sessInd}.L.ecc;
eccTransR = fData.voxProp{subjInd,sessInd}.R.eccTrans;
eccR = fData.voxProp{subjInd,sessInd}.R.ecc;
eccMax = max(fData.voxProp{subjInd,sessInd}.L.ecc);
Rmax = max(fData.voxProp{subjInd,sessInd}.R.ecc);

tickLabel = 0:0.2:eccMax;
tickLabel2 = 0:1:eccMax;
ax(2).YTick = eccTrans.toFlat{3}(tickLabel)';
ax(2).YTickLabel = cellstr(num2str(tickLabel'));
ax(2).YTickLabel(~ismembertol(tickLabel,tickLabel2)) = {''};
ax(2).YTickLabel(ismembertol(tickLabel,tickLabel2)) = cellstr(num2str(tickLabel2'));
ax(2).TickDir = 'out';

tickLabel = 0:0.2:Rmax;
tickLabel2 = 0:1:Rmax;
ax(1).YTick = flip(-eccTransR.toFlat{3}(tickLabel)',2);
ax(1).YTickLabel = flip(cellstr(num2str(tickLabel')),1);
ax(1).YTickLabel(flip(~ismembertol(tickLabel,tickLabel2),2)) = {''};
ax(1).YTickLabel(flip(ismembertol(tickLabel,tickLabel2),2)) = flip(cellstr(num2str(tickLabel2')),1);
ax(1).TickDir = 'out';

ax(1).XAxis.Visible = 'off';
ax(2).XAxis.Visible = 'off';

% Colorbar
hCb = colorbar;
hCb.Ticks = linspace(0,6,7)/6.*pi;
hCb.TickLabels = cellstr(num2str(linspace(0,6,7)'));
hCb.TickDirection = 'out';


% Voxel density
tmp2 = load(fullfile(funPath,inDir,[p.meta.subjList{subjInd} '.mat']));
res = tmp2.res; clear tmp2
d = res.(['sess' num2str(sessInd)]);

[eccLbefore,eccRbefore,densityLbefore,densityRbefore] = compareEcc(d);
[eccLafter,eccRafter,densityLafter,densityRafter] = compareEcc(d,eccTrans,eccTransR);

fDensity = figure('WindowStyle','docked');
tiledlayout(1,2,'TileIndexing','columnmajor')
nexttile

indLbefore = ~isnan(densityLbefore);
densityLbefore(indLbefore) = densityLbefore(indLbefore).*mode(diff(eccLbefore(indLbefore)));
indLafter = ~isnan(densityLafter);
densityLafter(indLafter) = densityLafter(indLafter).*mode(diff(eccLafter(indLafter)));

indRbefore = ~isnan(densityRbefore);
densityRbefore(indRbefore) = densityRbefore(indRbefore).*mode(diff(eccRbefore(indRbefore)));
indRafter = ~isnan(densityRafter);
densityRafter(indRafter) = densityRafter(indRafter).*mode(diff(eccRafter(indRafter)));


hAreaLbefore = area(-eccLbefore,densityLbefore./nansum(densityLbefore)); hold on
hAreaLafter = area(-eccLafter,densityLafter./nansum(densityLafter));
axL = gca;
axL.XTick = -flip(p.featSel.fov.threshVal);
axL.YAxis.Visible = 'off';
nexttile
hAreaRbefore = area(eccRbefore,densityRbefore./nansum(densityRbefore)); hold on
hAreaRafter = area(eccRafter,densityRafter./nansum(densityRafter));
axR = gca;
axR.XTick = p.featSel.fov.threshVal;
axR.YAxis.Visible = 'off';

hAreaLbefore.FaceAlpha = 0.5;
hAreaLafter.FaceAlpha = 0.5;

hAreaRbefore.FaceAlpha = 0.5;
hAreaRafter.FaceAlpha = 0.5;
