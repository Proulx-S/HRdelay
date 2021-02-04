function res = runAllDecoding(figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end

if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end

svmSpace = 'cart';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartNoAmp';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartNoDelay';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

svmSpace = 'cartReal';
dataType = 'wave';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);
dataType = 'sin';
[resBS.([svmSpace '_' dataType]),resWS.([svmSpace '_' dataType])] = runDecoding(svmSpace,dataType,verbose);

% [P,~,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right');
% disp('cart VS cartReal:')
% disp(['signed rank = ' num2str(STATS.signedrank)])
% disp(['one-sided p = ' num2str(P)])

disp('Between-session')
plotDecoding_acc(resBS,figOption,verbose)
disp('Witheen-session')
plotDecoding_acc(resWS,figOption,verbose)

res.BS = resBS;
res.WS = resWS;


% resBS.cartNoAmp.sess.acc
% 1-resBS.cartNoAmp.sess.acc_p