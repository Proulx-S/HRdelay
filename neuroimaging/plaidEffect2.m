clear all
close all
% cd('C:\Users\Sebastien\OneDrive - McGill University\work\projects\170210_HRdecoding\C_processing')
cd('C:\Users\sebas\OneDrive - McGill University\McGill\work\projects\170210_HRdecoding\C_processing')


fileList1 = {'02jp_sess1' '03sk_sess1' '04sp_sess1' '05bm_sess1' '06sb_sess1' '07bj_sess1'};
fileList2 = {'02jp_sess2' '03sk_sess2' '04sp_sess2' '05bm_sess2' '06sb_sess2' '07bj_sess2'};

subjInd = 3;
sessInd = 1;

for i = 1:length(fileList1)
    fullfile(pwd,[fileList1{i} '.mat'])
    load(fullfile(pwd,[fileList1{i} '.mat']),'d','resp')
    ori1{i,1} = d.xDataFixed(1:end/2,:);
    ori1r{i,1} = d.xData(1:end/2,:);
    ori2{i,1} = d.xDataFixed(end/2+1:end,:);
    ori2r{i,1} = d.xData(end/2+1:end,:);
    plaid{i,1} = d.normDataFixed;
    
    respOri1(:,i,1) = mean(resp.data.fix(:,:,1),1);
    respOri2(:,i,1) = mean(resp.data.fix(:,:,2),1);
    respPlaid(:,i,1) = mean(resp.data.fix(:,:,3),1);
    
    
    

    
    load(fullfile(pwd,[fileList2{i} '.mat']),'d','resp')
    ori1{i,2} = d.xDataFixed(1:end/2,:);
    ori2{i,2} = d.xDataFixed(end/2+1:end,:);
    plaid{i,2} = d.normDataFixed;
    
    respOri1(:,i,2) = mean(resp.data.fix(:,:,1),1);
    respOri2(:,i,2) = mean(resp.data.fix(:,:,2),1);
    respPlaid(:,i,2) = mean(resp.data.fix(:,:,3),1);
end

respOri1 = mean(respOri1,3).*100;
respOri2 = mean(respOri2,3).*100;
respPlaid = mean(respPlaid,3).*100;
%shift response
respOri1 = circshift(respOri1,-3,1);
respOri2 = circshift(respOri2,-3,1);
respPlaid = circshift(respPlaid,-3,1);



%%%%
%% Cartesian mean for orientations; Cartesian mean for voxels and Cartesian mean for sessions !!!!
%%%%
for subj = 1:size(ori1,1)
    ori_cart(subj,1) = mean(mean(cat(3,mean(ori1{subj,1},1),mean(ori2{subj,1},1)),3),2);
    ori_cart(subj,2) = mean(mean(cat(3,mean(ori1{subj,2},1),mean(ori2{subj,2},1)),3),2);
    
    plaid_cart(subj,1) = mean(mean(plaid{subj,1},1),2);
    plaid_cart(subj,2) = mean(mean(plaid{subj,2},1),2);
end
ori_cart = mean(ori_cart,2);
plaid_cart = mean(plaid_cart,2);

%mag
ori_mag = abs(ori_cart);
plaid_mag = abs(plaid_cart);

figure('WindowStyle','docked');
bar(pseudomedian([plaid_mag ori_mag], 1),'edgeColor','k','faceColor','w','lineWidth',5); hold on
plot([plaid_mag ori_mag]','o');
set(gca,'xTick',[1 2]); set(gca,'xTickLabel',{'plaid' 'ori'});
ylabel('amp')
[P_nonParam,H,STATS] = signrank(plaid_mag,ori_mag);
[H,P_param,CI,STATS] = ttest(plaid_mag,ori_mag);
title(['delta=' num2str(diff(pseudomedian([plaid_mag ori_mag], 1))) '; p(nonParam)=' num2str(P_nonParam) '; p(param)=' num2str(P_param)]);

%phase
ori_phase = angle(ori_cart);
plaid_phase = angle(plaid_cart);

figure('WindowStyle','docked');
bar(pseudomedian(wrapToPi([plaid_phase ori_phase]*-1+pi/2), 1)./pi*6,'edgeColor','k','faceColor','w','lineWidth',5); hold on
plot(wrapToPi([plaid_phase ori_phase]'*-1+pi/2)./pi*6,'o'); %ylim([5 8])
set(gca,'xTick',[1 2]); set(gca,'xTickLabel',{'plaid' 'ori'});
ylabel('delay (sec)')
[P_nonParam,H,STATS] = signrank(plaid_phase,ori_phase);
[P_param, F] = circ_htest(plaid_phase,ori_phase); % parametric; for paired samples
title(['delta=' num2str(diff(pseudomedian([plaid_phase ori_phase], 1)./pi*6 *-1+6-3)) '; p(nonParam)=' num2str(P_nonParam) '; p(param)=' num2str(P_param)]);

figure('WindowStyle','docked');
plot(wrapToPi([plaid_phase ori_phase]'*-1+pi/2)./pi*6); %ylim([5 8])




%% Plot responses
respOri = mean(cat(3,respOri1,respOri2),3);
figure('WindowStyle','docked');
errorbar(0:11,mean(respOri,2),std(respOri,[],2)./sqrt(size(respOri,2)),'r'); hold on
errorbar(0:11,mean(respPlaid,2),std(respPlaid,[],2)./sqrt(size(respPlaid,2)),'k');
legend(char({'Mean of responses to both stimulus orientations' 'Response to the mean of both stimulus orientations'}))


% %normalize according to mag
% respOri2 = respOri./repmat(mean([plaid_mag ori_mag],2)',[12 1]).*repmat(mean(mean([plaid_mag ori_mag],2),1),[12 6]);
% respPlaid2 = respPlaid./repmat(mean([plaid_mag ori_mag],2)',[12 1]).*repmat(mean(mean([plaid_mag ori_mag],2),1),[12 6]);
% figure('WindowStyle','docked');
% errorbar(0:11,mean(respOri2,2),std(respOri2,[],2)./sqrt(size(respOri2,2)),'r'); hold on
% errorbar(0:11,mean(respPlaid2,2),std(respPlaid2,[],2)./sqrt(size(respPlaid2,2)),'k');
% legend({'Mean of responses to both stimulus orientations' 'Response to the mean of both stimulus orientations'})


figure('WindowStyle','docked');
errorbar(0:11,mean(respOri1,2),std(respOri1,[],2)./sqrt(size(respOri1,2)),'r'); hold on
errorbar(0:11,mean(respOri2,2),std(respOri2,[],2)./sqrt(size(respOri2,2)),'b');
errorbar(0:11,mean(respPlaid,2),std(respPlaid,[],2)./sqrt(size(respPlaid,2)),'k');
legend(char({'ori1' 'ori2' 'plaid'}))

tmp = mean(respOri,2)
tmp = std(respOri,[],2)./sqrt(size(respOri,2))
tmp = mean(respPlaid,2)
tmp = std(respPlaid,[],2)./sqrt(size(respPlaid,2))


% [pval, F] = circ_htest(plaid_phase,ori_phase); % parametric; for paired samples
% [pval, table] = circ_wwtest(plaid_phase, ori_phase); % parametric; does not seem to be a test for paired samples
% pval = circ_cmtest(plaid_phase, ori_phase); % parametric; definitively not a test for paired samples
% circ_mtest
% [pval, z] = circ_mtest(plaid_phase-ori_phase, 0)
% [pval] = circ_medtest(plaid_phase-ori_phase,0); % boils down to a binomial test



for vox = 1:size(ori1r{subjInd,sessInd},2)
    [pval(vox), ~] = circ_wwtest(ori1r{subjInd,sessInd}(:,vox),ori2r{subjInd,sessInd}(:,vox));
end
[a,vox] = min(pval)

tmp1 = [abs(ori1r{subjInd,sessInd}(:,vox)) (angle(ori1r{subjInd,sessInd}(:,vox))./pi*6 *-1 +6-3)]
tmp2 = [abs(ori2r{subjInd,sessInd}(:,vox)) (angle(ori2r{subjInd,sessInd}(:,vox))./pi*6 *-1 +6-3)]
tmp1 = [real(ori1r{subjInd,sessInd}(:,vox)) imag(ori1r{subjInd,sessInd}(:,vox))]
tmp2 = [real(ori2r{subjInd,sessInd}(:,vox)) imag(ori2r{subjInd,sessInd}(:,vox))]

% [tmpX,tmpY] = pol2cart(angle(ori1{subjInd,sessInd})+pi/2,abs(ori1{subjInd,sessInd}));
clear tmpX1 tmpY1 tmpX2 tmpY2
[tmpX1,tmpY1] = pol2cart(wrapToPi(angle(plaid{subjInd,1})*-1+pi/2),abs(plaid{subjInd,1}));
[tmpX2,tmpY2] = pol2cart(wrapToPi(angle(plaid{subjInd,2})*-1+pi/2),abs(plaid{subjInd,2}));
angle(mean([mean(complex(tmpX1,tmpY1)) mean(complex(tmpX2,tmpY2))]))/pi*6

[a,b] = sort(exampleF,'descend')
ind = b; %(1:end/2)
% ind = randperm(size(tmpX1,2),round(size(tmpX1,2)/4))
angle(mean(complex(tmpX1(ind),tmpY1(ind))))/pi*6
figure('WindowStyle','docked');
gray = 0
h = scatter(tmpX1(ind),tmpY1(ind),ones(size(tmpY1(ind))).*5,'o','markerFaceColor','k','markerEdgeColor','none'); hold on
alpha(h,0.4)
axis equal
r=1;
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp,yp,'color',[gray gray gray]);
lims = [-max([abs(tmpX1(ind)) abs(tmpY1(ind))]) max([abs(tmpX1(ind)) abs(tmpY1(ind))])];
xlim(lims./3)
ylim(lims./3)
plot(get(gca,'xlim'),[0 0],'color',[gray gray gray])
plot([0 0],get(gca,'ylim'),'color',[gray gray gray])
plot([0 mean(tmpX1(ind))]*10,[0 mean(tmpY1(ind))]*10,'r')
plot([0 -mean(tmpX1(ind))]*10,[0 -mean(tmpY1(ind))]*10,'r')
set(gca, 'Color', 'None')
set(gcf, 'Color', 'White')
set(gca,'XTick',[]); set(gca,'XColor','none')
set(gca,'YTick',[]); set(gca,'YColor','none')
print('ScreenSizeFigure','-dpng','-r0')


% angle(ori1{subjInd,1})

% magRange = 0.2;
% magInd = abs(ori1{subjInd,1}) > abs(mean(ori1{subjInd,1}))-magRange  &  abs(ori1{subjInd,1}) < abs(mean(ori1{subjInd,1}))+magRange;
% find(magInd)
% phaseRange = 0.4;
% phaseInd = angle(ori1{subjInd,1})/pi*6 > angle(mean(ori1{subjInd,1}))/pi*6-phaseRange & angle(ori1{subjInd,1})/pi*6 < angle(mean(ori1{subjInd,1}))/pi*6+phaseRange;
% find(phaseInd)
% max(exampleF(magInd & phaseInd))

%get a descent voxel
Find = exampleF>40 & exampleF<60;
magTmp = abs(ori1{subjInd,1}(Find));
phaseTmp = angle(ori1{subjInd,1}(Find))/pi*6;
fTmp = exampleF(Find);
figure('WindowStyle','docked');
scatter3(magTmp,phaseTmp,fTmp); hold on
xlabel('mag'); ylabel('phase'); zlabel('F')
vox = find(cursor_info.Position(1)==magTmp)
find(cursor_info.Position(2)==phaseTmp)
tmp = find(Find); vox = tmp(vox);
vox = 1417;
[abs(ori1{subjInd,1}(vox)) abs(mean(ori1{subjInd,1}))]
[angle(ori1{subjInd,1}(vox))/pi*6 angle(mean(ori1{subjInd,1}))/pi*6]
scatter3(abs(ori1{subjInd,1}(vox)),angle(ori1{subjInd,1}(vox))/pi*6,exampleF(vox),'rx');

%plot time course
tmp = mean(exampleResp(vox,:,:,1),3);
t = (0:95)/6*pi;
[fitresult, gof] = sinFit(t, tmp, 0);
figure('WindowStyle','docked');
tmpX = (0:95)'; tmpY = tmp';
plot(tmpX,tmpY,'k'); hold on
tmpX = linspace(0,95,1000)'; tmpY = fitresult(linspace(t(1),t(end),1000));
plot(tmpX,tmpY,'r');
tmpY = tmp'*100;
tmpY = fitresult(t)*100;

size(exampleResp.data)
size(resp.data)

ind(1)

