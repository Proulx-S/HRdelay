function [resBS,resWS] = runAllDecoding(p,figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end

if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
clear resBS resWS
p.dataType = 'sin';
p.svmSpace = 'cartReal';
p.svmSpace = 'cartReal_affineRot';
p.svmSpace = 'cartImag';
p.svmSpace = 'cartImag_affineRot';
p.svmSpace = 'cart';
p.svmSpace = 'cartNoAmp';
p.svmSpace = 'cartNoAmpImag';
p.svmSpace = 'cartNoAmpImag_affineRot';
p.svmSpace = 'cartNoDelay';
p.svmSpace = 'cart_roi';
p.svmSpace = 'cart_affineRot';
p.svmSpace = 'cartNoAmp_affineRot';
p.svmSpace = 'cartNoAmp_affineRot_affineCart';
p.condPair = 'grat1VSgrat2';


p.condPair = 'grat1VSplaid';
% p.condPair = 'grat2VSplaid';
% p.condPair = 'grat1VSgrat2';
p.svmSpace = 'cartNoDelay';
[resBS,resWS] = runDecoding(p,verbose);
return

p.dataType = 'sin';
p.svmSpace = 'cart';
p.svmSpace = 'cartNoAmp';
p.svmSpace = 'cartNoDelay';
p.svmSpace = 'cartReal';
p.condPair = 'grat1VSgrat2';
p.condPair = 'grat1VSplaid';
p.condPair = 'grat2VSplaid';
[resBS.x,resWS.x] = runDecoding(p,verbose);

[resBS.polRoi_cartRoi,resWS.polRoi_cartRoi] = runDecoding(p,verbose);
[resBS.polRoi_cartVox,resWS.polRoi_cartVox] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoVoxThetaRoi_cartVox,resWS.polRhoVoxThetaRoi_cartVox] = runDecoding(svmSpace,dataType,verbose);
[resBS.polVox_cartVox,resWS.polVox_cartVox] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoRoiThetaVox_cartVox,resWS.polRhoRoiThetaVox_cartVox] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoVoxThetaVox_cartNone,resWS.polRhoVoxThetaVox_cartNone] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoVoxThetaVox_cartShiftRoiScaleNone,resWS.polRhoVoxThetaVox_cartShiftRoiScaleNone] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoVoxThetaVox_cartShiftRoiScaleRoi,resWS.polRhoVoxThetaVox_cartShiftRoiScaleRoi] = runDecoding(svmSpace,dataType,verbose);
[resBS.polRhoVoxThetaVox_cartShiftVoxScaleVox,resWS.polRhoVoxThetaVox_cartShiftVoxScaleVox] = runDecoding(svmSpace,dataType,verbose);

figure('WindowStyle','docked');
compareRes(resBS.polRhoRoiThetaVox_cartVox,resBS.polRhoVoxThetaRoi_cartVox)
figure('WindowStyle','docked');
compareRes(resBS.polRoi_cartVox,resBS.polRhoVoxThetaRoi_cartVox)

figure('WindowStyle','docked');
subplot(2,1,1)
compareRes(resBS.polRhoVoxThetaRoi_cartVox,resBS.polVox_cartVox)
subplot(2,1,2)
compareRes(resWS.polRhoVoxThetaRoi_cartVox,resWS.polVox_cartVox)

figure('WindowStyle','docked');
subplot(2,1,1)
compareRes(resBS.polRoi_cartVox,resBS.polRhoRoiThetaVox_cartVox)
subplot(2,1,2)
compareRes(resWS.polRoi_cartVox,resWS.polRhoRoiThetaVox_cartVox)


dataType = 'wave';
svmSpace = 'cart';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoAmp';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoDelay';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartReal';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);





clear resBS resWS
dataType = 'waveFull';
svmSpace = 'cart';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoAmp';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoAmpImag';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoDelay';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartReal';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartImag';
[resBS,resWS] = runDecoding(svmSpace,dataType,verbose);


% dataType = 'waveFull';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave',dataType,svmSpace,'between-session')
% compareDataType(resWS,'wave',dataType,svmSpace,'within-session')
% dataType = 'waveRun';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave',dataType,svmSpace,'between-session')
% compareDataType(resWS,'wave',dataType,svmSpace,'within-session')
% dataType = 'waveTrialSparse';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave',dataType,svmSpace,'between-session')
% compareDataType(resWS,'wave',dataType,svmSpace,'within-session')





svmSpace = 'cart';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseCat2';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseRep';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartNoAmp';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseCat2';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseRep';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartNoDelay';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseCat2';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseRep';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartReal';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseCat2';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'waveTrialSparseRep';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);


svmSpace = 'cart';
compareDataType(resBS,'waveTrialSparseCat2','waveTrialSparseRep',svmSpace,'between-session')
svmSpace = 'cartNoAmp';
compareDataType(resBS,'waveTrialSparseCat2','waveTrialSparseRep',svmSpace,'between-session')
svmSpace = 'cartNoDelay';
compareDataType(resBS,'waveTrialSparseCat2','waveTrialSparseRep',svmSpace,'between-session')
svmSpace = 'cartReal';
compareDataType(resBS,'waveTrialSparseCat2','waveTrialSparseRep',svmSpace,'between-session')


close all
clear resBS resWS
dataType = 'wave';
svmSpace = 'cart';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoAmp';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartNoDelay';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
svmSpace = 'cartReal';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% not scaled at svm; roi-wise rho scale
aBS = resBS;
aWS = resWS;
% scaled at svm; roi-wise rho scale
bBS = resBS;
bWS = resWS;
% not scaled at svm; voxel-spec rho scale
cBS = resBS;
cWS = resWS;
% scaled at svm; voxel-spec rho scale
dBS = resBS;
dWS = resWS;
save tmp aBS aWS bBS bWS cBS cWS dBS dWS

plotDecoding_acc(aBS,figOption,verbose)
plotDecoding_acc(bBS,figOption,verbose)

figure('WindowStyle','docked');
compareRes(aBS.cart_wave,bBS.cart_wave)
figure('WindowStyle','docked');
compareRes(bBS.cart_wave,dBS.cart_wave)
figure('WindowStyle','docked');
compareRes(aBS.cart_wave,cBS.cart_wave)

% svmSpace = 'cart';
% dataType = 'wave';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% dataType = 'sin';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave','sin',svmSpace,'between-session')
% compareDataType(resWS,'wave','sin',svmSpace,'within-session')
% 
% 
% svmSpace = 'cartNoAmp';
% dataType = 'wave';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% dataType = 'sin';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave','sin',svmSpace,'between-session')
% compareDataType(resWS,'wave','sin',svmSpace,'within-session')
% 
% svmSpace = 'cartNoDelay';
% dataType = 'wave';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% dataType = 'sin';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave','sin',svmSpace,'between-session')
% compareDataType(resWS,'wave','sin',svmSpace,'within-session')
% 
% svmSpace = 'cartReal';
% dataType = 'wave';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% dataType = 'sin';
% [resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
% compareDataType(resBS,'wave','sin',svmSpace,'between-session')
% compareDataType(resWS,'wave','sin',svmSpace,'within-session')

% [P,~,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right');
% disp('cart VS cartReal:')
% disp(['signed rank = ' num2str(STATS.signedrank)])
% disp(['one-sided p = ' num2str(P)])

disp('Between-session')
tmpField = fields(resBS);
tmp = struct2cell(resBS);
tmp = cell2struct(tmp(1:2:end),tmpField(1:2:end));
plotDecoding_acc(tmp,figOption,verbose)
disp('Witheen-session')
tmpField = fields(resWS);
tmp = struct2cell(resWS);
tmp = cell2struct(tmp(1:2:end),tmpField(1:2:end));
plotDecoding_acc(tmp,figOption,verbose)

res.BS = resBS;
res.WS = resWS;


% resBS.cartNoAmp.sess.acc
% 1-resBS.cartNoAmp.sess.acc_p