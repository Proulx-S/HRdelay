clear all
%% Dependencies
addpath(genpath(fullfile(pwd,'matlabFun')));
verbose = 0; % prints more info

%% Between-session feature selection parameters
% Activated voxels
featSel_bSess.activation.doIt = 1;
featSel_bSess.activation.fitType = 'fixed';
featSel_bSess.activation.threshType = 'p';
featSel_bSess.activation.threshVal = 0.05;
% Vein voxels
featSel_bSess.vein.doIt = 1;
featSel_bSess.vein.source = 'fullModelResid';% 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
featSel_bSess.vein.percentile = 20;
% Discriminant voxels
featSel_bSess.discrim.doIt = 1;
%% Display parameters
figOption.save = 0; % save all figures
figOption.subj = 1; % subjInd-> plots participants subjInd; +inf-> plots all participant (if verbose==0, will only plot subjInd==1 but still produce and save all the other figures)

runFit(verbose)

preprocAndShowMasks(featSel_bSess,figOption,verbose)

inspectSubjAndExclude(figOption,verbose)

runGroupAnalysis_sin(figOption,verbose)

runGroupAnalysis_hr(figOption,verbose)

res = runAllDecoding(figOption,verbose);
nObs = sum(res.cart.nObs(:));
accCart = res.cart.summary.acc;
cart = zeros(1,nObs);
cart(round(1:accCart*nObs)) = 1;
accCartReal = res.cartReal.summary.acc;
cartReal = zeros(1,nObs);
cartReal(round(1:accCartReal*nObs)) = 1;
X = cart;
Y = cartReal;
[z,pv1,pv2,pv3]=compare_bino_prob(X,Y)
[z,pv1,pv2,pv3]=compare_bino_prob(Y,X)

hitCart = round(accCart*nObs);
hitCartReal = round(accCartReal*nObs);
x = [hitCart hitCartReal
    nObs-hitCart nObs-hitCartReal];
STATS=mybarnard(x,'Tbx',1000)

prob = binopdf(0:nObs,nObs,accCart)'*binopdf(0:nObs,nObs,accCartReal);
imagesc(prob)
ax = gca;
ax.PlotBoxAspectRatio = [1 1 1];
ax.YDir = 'normal';
hold on
plot([0 nObs],[0 nObs],'r')
plot([nObs 0],[0 nObs],'r')

cartLarger = 0;
cartSmaller = 0;
cartEqual = 0;
for i = 1:nObs
    for ii = 1:nObs
        if i>ii
            cartLarger = cartLarger + prob(i,ii);
        elseif i<ii
            cartSmaller = cartSmaller + prob(i,ii);
        else
            cartEqual = cartEqual + prob(i,ii);
        end
    end
end
1-(cartLarger)


return
runAllDecodingPerm(res,figOption,verbose);
