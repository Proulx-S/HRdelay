function f = plotAllDecoding3(p,res,info)
f1 = plotAllDecoding2(p,res,info);
f2 = plotAllDecoding(p,res,info);

delta = 2;
ax1 = findobj(f1.Children,'Type','Axes');
ax2 = findobj(f2.Children,'Type','Axes');

hBar1 = findobj(ax1.Children,'Type','Bar');
hBar2 = findobj(ax2.Children,'Type','Bar');


figure(f1)
h = hBar1(1);
D_x = mode(diff(h.XData));
wTx = D_x;
b_w = h.BarWidth;
bg_w = b_w;
b_x = b_w*wTx;
bg_x = bg_w*wTx;

wTx_1 = wTx;
D_x1 = D_x;
b_x1 = b_x;
bg_x1 = bg_x;


h = hBar2(1);
D_x = mode(diff(h.XData));
wTx = D_x/3.5;
b_w = h.BarWidth;
bg_w = b_w*2+(1-b_w);
b_x = b_w*wTx;
bg_x = bg_w*wTx;

wTx_2 = wTx;
D_x2 = D_x;
b_x2 = b_x;
bg_x2 = bg_x;


% figure(f1)
D1 = D_x1;
b1 = b_x1;
bg1 = bg_x1;
d1 = D1-bg1;
% figure(f2)
D2 = D_x2;
b2 = b_x2;
bg2 = bg_x2;
d2 = D2-bg2;

%D1=d2+b2
scale_1 = (d2+b2)/D1;
wTx_1 = wTx_1*scale_1;
d1 = d1*scale_1;
%b1=bg2
b1_w = b2/wTx_1;

figure(f2)
hTmp = copyobj(findobj(ax1.Children,'Type','Bar'),ax2);
hTmp.XData = hTmp.XData*scale_1 - delta;
hTmp.BarWidth = b1_w;
hBar1 = hTmp;
hBar2 = findobj(ax2.Children,'Type','Bar');

hTmp = copyobj(findobj(ax1.Children,'Type','ErrorBar'),ax2);
hTmp.XData = (hTmp.XData*scale_1 - delta);
hEr1 = hTmp;
hEr2 = findobj(ax2.Children,'Type','ErrorBar');

hTmp = copyobj(findobj(ax1.Children,'Type','Patch'),ax2);
for i = 1:length(hTmp)
    hTmp(i).Vertices(:,1) = hTmp(i).Vertices(:,1)*scale_1 - delta + (b1-b2)/4;
end
hPatch1 = hTmp;
hPatch2 = findobj(ax2.Children,'Type','Patch');

hTmp = copyobj(findobj(ax1.Children,'Type','Scatter'),ax2);
for i = 1:length(hTmp)
    hTmp(i).XData = hTmp(i).XData*scale_1 - delta;
end
hScat1 = hTmp;
hScat2 = findobj(ax2.Children,'Type','Scatter');



fac = 0.5;
D1 = mode(diff(hBar1.XData));
% hBar1.XData
delta = (2:-1:0)*(d1*fac/D1);
% fac2 = d1*fac/D1;
hBar1.XData = hBar1.XData + delta;
fac2 = D1/mode(diff(hBar1.XData));
hBar1.BarWidth = hBar1.BarWidth*fac2;
hEr1.XData = hEr1.XData + delta;
for i = 1:length(hScat1)
    hScat1(i).XData = hScat1(i).XData + delta;
end
for i = 1:length(hPatch1)
    hPatch1(i).Vertices(:,1) = hPatch1(i).Vertices(:,1) + delta(i);
end
uistack(hPatch1,'bottom')
% set(hPatch1,'FaceColor','r')


delete(f1)
f = f2; clear f2
ax = findobj(f.Children,'Type','Axes');
hLine = findobj(ax.Children,'Type','Line');
hLine.XData = xlim;
uistack(hLine,'bottom');

hScat = findobj(ax.Children,'Type','Scatter');
set(hScat,'MarkerEdgeColor',[1 1 1].*0.7);
tmp = get(hScat,'SizeData');
set(hScat,'SizeData',mode([tmp{:}])/2);

delete(findobj(f.Children,'Type','Legend'));
delete(findobj(f.Children,'Type','Text'));
ax.Box = 'off';
ax.Color = 'none';
f.Color = 'w';




%% Save
if p.figOption.save
    fullfilename = fullfile(p.figOption.finalDir,'decodingSummary3');
    curF = f;
    curF.Color = 'none';
    set(findobj(curF.Children,'type','Axes'),'color','none')
    curFile = fullfilename;
    curExt = 'eps';
    exportgraphics(curF,[curFile '.' curExt],'ContentType','vector'); if p.figOption.verbose; disp([curFile '.' curExt]); end
    curExt = 'svg';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    curF.Color = 'w';
    curExt = 'fig';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    curExt = 'jpg';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
end

function [yBoot_90CI,yBoot_95CI] = bootThat(y,nBoot)
nSubj = size(y,3);
sz = size(y);
y = permute(y,[length(sz) 1:(length(sz)-1)]);
sz = sz([length(sz) 1:(length(sz)-1)]);
sz(1) = nBoot;
yBoot = nan(sz);
for bootInd = 1:nBoot
    i = nan(nSubj,1);
    for drawInd = 1:nSubj
        i(drawInd) = randperm(nSubj,1);
    end
    yBoot(bootInd,:) = mean(y(i,:),1);
end

yBoot = permute(yBoot,[2:length(sz) 1]);
sz = sz([2:length(sz) 1]);
yBoot_90CI = prctile(yBoot,[5 95],length(sz));
yBoot_95CI = prctile(yBoot,[2.5 97.5],length(sz));


