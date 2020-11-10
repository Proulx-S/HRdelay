function validation
%% Define paths
repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataDir = 'C-derived\DecodingHR';
anatPath = fullfile(repoPath,dataDir,'anat');
anatLevel = 'z';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel1 = 'x';
funLevel2 = 'y';
funLevel3 = 'z';

fileList = {'02jp_sess1' '03sk_sess1' '04sp_sess1' '05bm_sess1' '06sb_sess1' '07bj_sess1';...
            '02jp_sess2' '03sk_sess2' '04sp_sess2' '05bm_sess2' '06sb_sess2' '07bj_sess2'}';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';


subjInd = 1;
%% Load fun data
load(fullfile(funPath,funLevel2,subjList{subjInd},'v1SinCos_1perRun_move12.mat'),'results')
% brain
tmpFilename = ls(fullfile(funPath,funLevel1,subjList{subjInd},'trun101_preprocessed.nii*'));
a = load_nii(fullfile(funPath,funLevel1,subjList{subjInd},tmpFilename));
brain = flipdim(permute(a.img,[3 1 2 4]),1); clear a
brain = brain(:,:,:,1);

%% Set masks
% maskV1
tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
maskV1 = flipdim(permute(a.img,[3 1 2 4]),1); clear a
maskV1(:,:,1) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
maskV1(:,:,end) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices

% maskECC
tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'lh.ecc.nii*'));
a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'rh.ecc.nii*'));
a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
eccR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
find((eccL~=0 & eccR~=0));
ecc = eccL; ecc(eccR~=0) = eccR(eccR~=0); clear eccL eccR
maskECC = ecc>1 & ecc<6;

% maskAnat
maskAnat = maskV1 & maskECC;

% maskFit area
maskFit = results.mask;

% maskFun 
for sessInd = 1:2
    eval(['F = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F(maskAnat(maskFit));']);
    eval(['df = results.OLS.mixed.Fsess' num2str(sessInd) '.df;']);
    
    [FDR,P] = getPfromF(F,df);
    tmp = nan(size(brain)); tmp(maskAnat) = F;
    F = tmp;
    
    tmp = nan(size(brain)); tmp(maskAnat) = FDR;
    FDR = tmp;
    
    tmp = nan(size(brain)); tmp(maskAnat) = P;
    P = tmp;
    
    maskFun.(['sess' num2str(sessInd)]).F = F;
    maskFun.(['sess' num2str(sessInd)]).P = P;
    maskFun.(['sess' num2str(sessInd)]).FDR = FDR;
end

%% Plot masking
figure();
ax1 = subplot(2,2,1);
imagesc(brain(:,:,10)); colormap gray; axis off;
ax1.PlotBoxAspectRatio = [1 1 1];
title('BOLD image (1 TR)')

ax2 = subplot(2,2,2);
imagesc(maskAnat(:,:,10)); colormap gray; axis off
ax2.PlotBoxAspectRatio = [1 1 1];
title({'Anatomical ROI' '\fontsize{8} V1, 1 to 6 dva eccentricity'})

ax3 = subplot(2,2,3); axis off
im1 = brain(:,:,10);
im1 = ind2rgb(uint8(im1./max(im1(:))*255),gray(256));
F = maskFun.sess1.F(:,:,10);
thresh = maskFun.sess1.P(:,:,10);
F(thresh>0.05) = 0; % functional threshold here
F(isnan(F)) = 0; % functional threshold here
im2 = ind2rgb(uint8(F./max(F(:))*255),autumn(256));
ind = repmat(F~=0,[1 1 3]);
im1(ind) = im2(ind);
imshow(im1)
ax3.PlotBoxAspectRatio = [1 1 1];

axCB = copyobj(ax3,gcf);
axes(axCB);
CB = permute(flipud(autumn(size(ax1.Children.CData,1))),[1 3 2]);
CB = repmat(CB,[1 round(0.1*size(ax1.Children.CData,2)) 1]);
imshow(CB);
axCB.Position(1) = ax3.Position(1) - ax3.Position(3)*0.5;
axCB.YAxis.Visible = 'on';
axCB.YAxis.TickValues = [1 size(ax1.Children.CData,1)];
tmp = 0.05-thresh(:); tmp(tmp<0) = nan; [~,b] = min(tmp);
axCB.YAxis.TickLabels = {num2str(max(F(:)),'%0.1f') num2str(F(b),'%0.1f')};
axCB.YAxis.Label.String = {'F' 'p<0.05'};
axCB.YAxis.Label.Rotation = 0;
axCB.YAxis.Label.VerticalAlignment = 'middle';
axCB.YAxis.Label.HorizontalAlignment = 'center';

ax3.Title.String = 'Visual Activation';

ax4 = subplot(2,2,4);
hm = histogram(maskFun.sess1.F(maskAnat),1000); hold on
xlim([0 20])
ylim(ylim*1.05)
hm.FaceColor = 'k';
ax4.PlotBoxAspectRatio = [1 1 1];
ylabel('vox count');
xlabel('F')
line([1 1].*F(b),ylim,'Color','r');
text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,'p=0.05')
title('All ROI voxels')
ax4.XTick = sort([ax4.XTick round(F(b),2)]);
ax4.TickDir = 'out';
ax4.XTick(ax4.XTick==0) = [];



%% Apply masks to fun data, vectorize voxels and split sessions and conditions
[X,Y] = pol2cart(results.OLS.mixed.delay,results.OLS.mixed.amp);
data = complex(X,Y); clear X Y

size(maskAnat(maskFit))
size(maskFit(maskAnat))

% mask and vectorize
sz = size(data);
dataX = nan(sz(4),sum(maskAnat(maskFit)));
for runInd = 1:sz(4)
    curData = data(:,:,:,runInd);
    tmp = false(size(curData));
    tmp(:) = maskAnat(maskFit);
    dataX(runInd,:) = curData(tmp);
end
list = {'F' 'FDR' 'P'};
for i = 1:3
    for sessInd = 1:2
        tmp = maskFun.(['sess' num2str(sessInd)]).(list{i});
        tmp = tmp(maskFit);
        tmp = tmp(maskAnat(maskFit))';
        stats.(list{i}).(['sess' num2str(sessInd)]) = tmp;
    end
end

% split sessions and conditions
sessLabel = [results.inputs.opt.sessionLabel{:}];
condLabel = repmat(1:3,[size(dataX,1)/3 1]); condLabel = condLabel(:)';
for sessInd = 1:2
    d.(['sess' num2str(sessInd)]).xData = nan(sum(sessLabel==sessInd)/3,size(dataX,2),3);
    for condInd = 1:3
        d.(['sess' num2str(sessInd)]).xData(:,:,condInd) = dataX(sessLabel==sessInd & condLabel==condInd,:);
    end
    for i = 1:3
        d.(['sess' num2str(sessInd)]).(list{i}) = stats.(list{i}).(['sess' num2str(sessInd)]);
    end
end








% new = cat(1,d.sess1.xData(:,:,1),d.sess1.xData(:,:,2),d.sess1.xData(:,:,3));
% old = load('C:\Users\sebas\OneDrive - McGill University\McGill\work\projects\170210_HRdecoding\C_processing\02jp_sess1.mat');
% old = cat(1,old.d.xData,old.d.normData);
% scatter(abs(old(:)),abs(new(:))); hold on
% line(xlim,ylim)



% 
% 
% 
% 
% 
% for sessInd = 1:2
%     eval(['Fsess' num2str(sessInd) ' = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F;']);
%     eval(['[FDRsess' num2str(sessInd) ',Psess' num2str(sessInd) '] = getPfromF(results.OLS.mixed.Fsess' num2str(sessInd) '.val.F,results.OLS.mixed.Fsess' num2str(sessInd) '.df);']);
% %     eval(['Fsess' num2str(sessInd) ' = single(zeros(size(brain)));']);
% %     eval(['Fsess' num2str(sessInd) '(maskFit) = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F;']);
% %     eval(['[tmpFDRsess' num2str(sessInd) ',tmpPsess' num2str(sessInd) '] = getPfromF(results.OLS.mixed.Fsess' num2str(sessInd) '.val.F,results.OLS.mixed.Fsess' num2str(sessInd) '.df);']);
% %     eval(['FDRsess' num2str(sessInd) ' = single(zeros(size(brain)));']);
% %     eval(['FDRsess' num2str(sessInd) '(maskFit) = tmpFDRsess' num2str(sessInd) ';']);
% %     eval(['Psess' num2str(sessInd) ' = single(zeros(size(brain)));']);
% %     eval(['Psess' num2str(sessInd) '(maskFit) = tmpPsess' num2str(sessInd) ';']);
% %     eval(['clear(''tmpFDRsess' num2str(sessInd) ''',''tmpPsess' num2str(sessInd) ''')'])
% end
% 
% % apply mask, vectorizing all data
% sz = size(data);
% maskFit = false(sz(1:3));
% maskFit(:) = maskAnat(maskFit);
% dataX = nan(size(data,4),sum(maskFit(:)));
% for i = 1:size(data,4)
%     tmp = data(:,:,:,i);
%     dataX(i,:) = tmp(maskFit);
% end
% for sessInd = 1:2
%     eval(['Fsess' num2str(sessInd) 'X = Fsess' num2str(sessInd) '(maskFit)'';']);
%     eval(['FDRsess' num2str(sessInd) 'X = FDRsess' num2str(sessInd) '(maskFit)'';']);
%     eval(['Psess' num2str(sessInd) 'X = Psess' num2str(sessInd) '(maskFit)'';']);
% end
% whos *sess1 *sess2 data
% whos *sess1X *sess2X dataX
% 
% % for sessInd = 1:2
% %     eval(['Fsess' num2str(sessInd) 'X = nan(1,sum(maskFit(:)));']);
% %     eval(['FDRsess' num2str(sessInd) 'X = nan(1,sum(maskFit(:)));']);
% %     eval(['Psess' num2str(sessInd) 'X = nan(1,sum(maskFit(:)));']);
% %     eval(['Fsess' num2str(sessInd) 'X(1,:) = Fsess' num2str(sessInd) '(maskFit);']);
% %     eval(['FDRsess' num2str(sessInd) 'X(1,:) = FDRsess' num2str(sessInd) '(maskFit);']);
% %     eval(['Psess' num2str(sessInd) 'X(1,:) = Psess' num2str(sessInd) '(maskFit);']);
% % end
% 
% 
% % Split sessions and cond
% sessList = [results.inputs.opt.sessionLabel{:}];
% condList = repmat(1:3,[size(dataX,1)/3 1]); condList = condList(:)';
% 
% for sessInd = 1:2
%     d.(['sess' num2str(sessInd)]).xData = nan(sum(sessList==sessInd)/3,size(dataX,2),3);
%     for condInd = 1:3
%         d.(['sess' num2str(sessInd)]).xData(:,:,condInd) = dataX(sessList==sessInd & condList==condInd,:);
%     end
%     d.(['sess' num2str(sessInd)]).F = eval(['Fsess' num2str(sessInd) 'X']);
%     d.(['sess' num2str(sessInd)]).FDR = eval(['FDRsess' num2str(sessInd) 'X']);
%     d.(['sess' num2str(sessInd)]).P = eval(['Psess' num2str(sessInd) 'X']);
%     d.(['sess' num2str(sessInd)]).info = '%BOLD: run x vox x cond[ori1,ori2,plaid]';
% end
% 
% 
% %% Plots
% % Show brain and masks
% % figure('WindowStyle','docked');
% % imagesc(brain(:,:,10)); colormap gray; axis square
% % % figure('WindowStyle','docked');
% % % imagesc(maskV1(:,:,10)); colormap gray; axis square
% % % figure('WindowStyle','docked');
% % % imagesc(ecc(:,:,10)); colormap gray; axis square
% % % figure('WindowStyle','docked');
% % % imagesc(maskECC(:,:,10)); colormap gray; axis square
% % figure('WindowStyle','docked');
% % imagesc(mask(:,:,10)); colormap gray; axis square
% % % figure('WindowStyle','docked');
% % % imagesc(maskFit(:,:,10)); colormap gray; axis square
% figure();
% ax1 = subplot(2,2,1);
% imagesc(brain(:,:,10)); colormap gray; axis square; axis off
% title('BOLD image (1 TR)')
% 
% ax2 = subplot(2,2,2);
% imagesc(maskAnat(:,:,10)); colormap gray; axis square; axis off
% title('V1 mask from 1 to 6 ecc')
% 
% ax3 = subplot(2,2,3);
% im1 = brain(:,:,10);
% im1 = ind2rgb(uint8(im1./max(im1(:))*255),gray(256));
% imTmp = nan(size(tmpMask));
% imTmp(tmpMask) = Fsess1(tmpMask);
% im2 = nan(size(brain(:,:,10)));
% im2(maskFit(:,:,10)) = imTmp;
% im2 = ind2rgb(uint8(im2./max(im2(:))*255),autumn(256));
% 
% 
% imshow(im2)
% imagesc(im2); colormap gray; axis square; axis off
% 
% im2 = ind2rgb(uint8(im1./max(im1(:))*255),gray(256));
% 
% imshow(im1)
% 
% im1 = ax1.Children.CData;
% im1 = uint8(im1./max(im1(:))*255);
% cmap1 = ax1.Colormap;
% im1 = ind2rgb(im1,cmap1);
% 
% tmpMask = maskFit(:,:,10);
% im2 = nan(size(tmpMask));
% im2(tmpMask) = Fsess1(tmpMask);
% im2(Psess1(:,:,10)>0.05) = nan;
% im2 = uint8(im2./max(im2(:)).*255);
% im2 = ind2rgb(im2,autumn(256));
% 
% ;
% figure('WindowStyle','docked');
% im2 = imagesc(im2); im2 = im2.CData;
% cmap2 = get(gca,'Colormap');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Fsess1(:,:,10)
% maskV1(:,:,10)
% im2(maskFit)
% whos Fsess1
% 
% 
% ax3 = subplot(2,2,3);
% imagesc(Fsess1(:,:,10),[]); axis square; axis off
% ax3 = gca;
% ax3.Colormap = colormap('autumn');
% 
% im2 = im1;
% hi = imagesc(Fsess1(:,:,10))
% hi
% % im2(maskFit(:,:,10)) = 
% imagesc(FDRsess1(:,:,10))
% 
% imshow(im1)
% 
% cmap1 = ax1.Children.C
% im1 = uint8(im1./max(im1(:))*255);
% im1 = ind2rgb(im1,gray(256));
% 
% Fsess1
% 
% imshow(im1)
% 
% 
% 
% Ftmp = nan(size(brain(:,:,10)));
% Ftmp(maskFit(:,:,10)) = Fsess1(:,:,10);
% FDRtmp = nan(size(brain(:,:,10)));
% FDRtmp(maskFit(:,:,10)) = FDRsess1(:,:,10);
% tmp = FDRtmp(:); [~,b] = min(abs(tmp-0.05));
% h = imagesc(Ftmp,[0 10]); axis square; axis off
% h.AlphaData = ~isnan(Ftmp);
% colormap(gca,'jet')
% hist(Ftmp(:))
% 
% imagesc(maskAnat(:,:,10)); axis square; axis off
% title('BOLD image (1 TR)')
% 
% subplot(2,2,4)
% imagesc(Fsess1(:,:,10)); axis square; axis off
% title('BOLD image (1 TR)')
% 
% subplot(2,2,3)
% imagesc(Fsess1(:,:,10)); axis square; axis off
% title('BOLD image (1 TR)')
% 
% 
% Fsess1
% imagesc(maskAnat(:,:,10))
% 
% % Show design matrices
% figure();
% imagesc(results.OLS.mixed.designmatrix(1:96,[1:2 67:80])'); colormap gray
% daspect([1 1 1])
% ylabel('regressors')
% xlabel('TRs')
% ax = gca; ax.XTick = []; ax.YTick = [];
% title('Design matrix for one run')
% 
% figure();
% nRun = sum(sessList==1);
% designmatrix = results.OLS.mixed.designmatrix(:,[1:(nRun*2) 66+(1:nRun*14)]);
% designmatrix(all(designmatrix==0,2),:) = [];
% imagesc(designmatrix'); colormap gray
% daspect([3 1 1])
% ylabel('regressors')
% xlabel('Runs and TRs')
% ax = gca; ax.XTick = []; ax.YTick = [];
% title('Full design matrix for one session')
% 
% 
% 
% 