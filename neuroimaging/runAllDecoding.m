function res = runAllDecoding(figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end

if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end

svmSpace = 'cart';
[resBS.(svmSpace),resWS.(svmSpace)] = runDecoding(svmSpace,verbose);

svmSpace = 'cartNoAmp';
[resBS.(svmSpace),resWS.(svmSpace)] = runDecoding(svmSpace,verbose);

svmSpace = 'polDelay';
[resBS.(svmSpace),resWS.(svmSpace)] = runDecoding(svmSpace,verbose);

svmSpace = 'cartNoDelay';
[resBS.(svmSpace),resWS.(svmSpace)] = runDecoding(svmSpace,verbose);

svmSpace = 'cartReal';
[resBS.(svmSpace),resWS.(svmSpace)] = runDecoding(svmSpace,verbose);

% [P,~,STATS] = signrank(mean(res.cart.acc,2),mean(res.cartReal.acc,2),'tail','right');
% disp('cart VS cartReal:')
% disp(['signed rank = ' num2str(STATS.signedrank)])
% disp(['one-sided p = ' num2str(P)])

disp('Between-session')
plotDecoding_acc(resBS,figOption,verbose)
disp('Between-session')
plotDecoding_acc(resWS,figOption,verbose)

res.BS = resBS;
res.WS = resWS;
